### Enable this line to reduce verbositys
# log.none()
from Bio.PDB import *
from .util import *
import sys,os

#### Comment these to recover warning on importing PDB
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)


_hydrogen = re.compile("[123 ]*H.*") 
class unsel_H(object): 

	"""	
	Adapted from Bio.PDB.Dice.ChainSelector(), internal routine for extract()
	Remove hydrogens, waters and ligands. Only use model 0 by default. 
	""" 
	
	def __init__(self, model_id = 0):
		self.model_id = model_id
	
	def accept_model(self, model): 
	    # model - only keep model 0 
	    if model.get_id() == self.model_id: 
	        return 1 
	    return 0 	
	def accept_chain(self, chain): 
	    return 1
	def accept_residue(self, residue): 
	    hetatm_flag, resseq, icode = residue.get_id() 
	    if hetatm_flag != " ": 
	        # skip HETATMS 
	        return 0 
	    if icode != " ": 
	        warnings.warn("WARNING: Icode %s at position %s" 
	                      % (icode, resseq), BiopythonWarning) 
	    return 1 
	
	def accept_atom(self, atom): 
	    # atoms - get rid of hydrogens 
	    name = atom.get_id() 
	    if _hydrogen.match(name): 
	        return 0 
	    else: 
	        return 1 
	
Hsel = unsel_H()
io = PDBIO()


def parse_PDB(pdbname,pdbdir=None,parser = None, **kwargs):
	if pdbdir:
		pass
	else:
		if 'PDBlib' in os.environ.keys():
			pdbdir = '$PDBlib/'
		else:
			os.environ['$PDBlib'] = './'
			pdbdir = '$PDBlib/'
	pdbfile = os.path.expandvars(pdbdir) + pdbname
	# print os.path.isfile(pdbfile)
	# assert os.path.isfile(pdbfile), "Cannot open PDB file at %s" % pdbfile
	
	if not parser:
		parser = PDBParser()
	struct = parser.get_structure( pdbname, pdbfile)

	return struct

if not os.path.isdir(full('$PDBlib/sanitised')):
	os.mkdir(full('$PDBlib/sanitised'))


def sanitise(struct, **kwargs):
	pdbname = struct.id
	san_pdbname = 'sanitised/%s'% pdbname
	san_pdbfile = full( '$PDBlib/' + san_pdbname)

	if os.path.isfile( san_pdbfile ):
		pass
	else: 
		struct = parse_PDB(pdbname, **kwargs)
		io.set_structure( struct )
		io.save(  san_pdbfile,
			Hsel)
	struct = parse_PDB(san_pdbname, **kwargs)
	return struct

def get_something(input, env = None, auto_complete = False, s0 = None, pdbdir = None, cutout=15.0,cutin=3.5,
	sanitising = 1,
	 **kwargs):
	if isinstance(input,Structure.Structure):
		struct = input
	else:
		pdbname = input
		struct = parse_PDB(pdbname, **kwargs)

	if sanitising:
		struct = sanitise(struct, **kwargs)

		# struct = p.get_structure('X', pdbfile)
	alst = list(struct.get_atoms())
	acount = len(alst)
	rcount = sum( 1 for _ in struct.get_residues())
	nbpair_count = (
	len( NeighborSearch(alst ).search_all(radius =cutout))  
	-  len( NeighborSearch(alst ).search_all(radius =cutin)) 
	)

	### For some weird reasons, many structures are not correctly parsed by Biopython, (but Okay with modeller)
	### This is a temporary patch to detect overlapping

	# if acount > rcount * 11:
	# 	acount /= 2
	# 	nbpair_count /= 4

	outdict = {"nDOPE":0,
		"DOPE": 0,
		"nbpair_count":nbpair_count,
		"atom_count": acount,
		 "res_count": rcount,
		  }
	return outdict


# if django.settings

if 'USE_MODELLER' in os.environ.keys():
	USE_MODELLER = int(os.environ['USE_MODELLER'])
else:
	USE_MODELLER = 0
if USE_MODELLER:
	from modeller import *
	from modeller.scripts import complete_pdb

	def	init_env(env=None):

		with stdoutIO() as s:	
			env = environ()
			#env.io.atom_files_directory = ['../atom_files']
			env.io.atom_files_directory = ['../pdbs','$(PDBlib)/',
			'$(repos)/cathdb/dompdb/',
			'$(repos)/cathdb/temppdbs/',
			]
			env.libs.topology.read(file='$(LIB)/top_heav.lib')
			env.libs.parameters.read(file='$(LIB)/par.lib')
		return env	

	def get_something_modeller( pdbfile, env = None, auto_complete = False, s0 = None, **kwargs):
		if not env:
			env = environ()
			#env.io.atom_files_directory = ['../atom_files']
			env.io.atom_files_directory = ['../pdbs',
			'$(PDBlib)']
			env.libs.topology.read(file='$(LIB)/top_heav.lib')
			env.libs.parameters.read(file='$(LIB)/par.lib')
		if auto_complete:
			mdl = complete_pdb(env, pdbfile)
		else:
			mdl = model(env)
			mdl.read( pdbfile, model_format='PDB', model_segment=('FIRST:@', 'LAST:'), io=None);
		
		if not s0:
			s0 = StringIO.StringIO();
		else:
			s0.truncate(0)
		
		with stdoutIO(s0) as s:
			nDOPE = mdl.assess_normalized_dope()
		s0buf = ''.join(s0.buflist)
		# return 
		outdict = {"nDOPE":nDOPE,
			"DOPE": float(p_energy.findall(s0buf)[0]),
			"nbpair_count":int(p_nb.findall(s0buf)[0]),
			"atom_count":int(p_atomCount.findall(s0buf)[0].strip()),
			 "res_count":  int(p_resCount.findall(s0buf)[0]),
			  }
		return outdict
else:
	def get_something_modeller( pdbfile, env = None, auto_complete = False, s0 = None, **kwargs):
		assert False, "'get_something_modeller()' is not defined!"

def get_nDOPE( pdbfile, env = None, auto_complete = False, **kwargs):
# pdbfile = "4xz8A_chop"
#mdl = complete_pdb(env, "1fas")
	
	### Use existing environment to avoid redundant re-initialisation.
	# if not env:
	# 	env = environ()
	# 	#env.io.atom_files_directory = ['../atom_files']
	# 	env.io.atom_files_directory = ['../pdbs']
	# 	env.libs.topology.read(file='$(LIB)/top_heav.lib')
	# 	env.libs.parameters.read(file='$(LIB)/par.lib')

	# if auto_complete:
	# 	mdl = complete_pdb(env, pdbfile)
	# else:
	# 	mdl = model(env)
	# 	mdl.read( pdbfile, model_format='PDB', model_segment=('FIRST:@', 'LAST:'), io=None);
	# nDOPE = mdl.assess_normalized_dope()
	# # pdbname = os.path.basename(pdbfile)


	nDOPE = get_something(pdbfile, env, auto_complete,**kwargs).get("nDOPE")
	return nDOPE