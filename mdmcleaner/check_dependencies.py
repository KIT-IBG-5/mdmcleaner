#!/usr/bin/env python
import subprocess
import sys
import re
import importlib

external_dependency_version_dict = {	"wget" : { "vpattern" : "Wget (\d+(\.\d+)+).*", "executable" : ["wget"], "version_option" : ["--version"], "required_version" : "1.19"},\
										"blastn": { "vpattern" : ": (\d+(\.\d+)+).*", "executable" : ["blastn"], "version_option" : ["-version"], "required_version" : "2.10"},\
										"blastp": { "vpattern" : ": (\d+(\.\d+)+).*", "executable" : ["blastp"], "version_option" : ["-version"], "required_version" : "2.10"},\
										"makeblastdb": { "vpattern" : ": (\d+(\.\d+)+).*", "executable" : ["makeblastdb"], "version_option" : ["-version"], "required_version" : "2.10"},\
										"blastdbcmd": { "vpattern" : ": (\d+(\.\d+)+).*", "executable" : ["blastdbcmd"], "version_option" : ["-version"], "required_version" : "2.10"},\
										"diamond": { "vpattern" : "diamond version (\d+(\.\d+)+).*", "executable" : ["diamond"], "version_option" : ["version"], "required_version" : "2.0.6"},\
										"hmmsearch": { "vpattern" : "# HMMER (\d+(\.\d+)+).*", "executable" : ["hmmsearch"], "version_option" : ["-h"], "required_version" : "3.3.1"},\
										"barrnap": { "vpattern" : "barrnap (\d+(\.\d+)+).*", "executable" : ["barrnap"], "version_option" : ["--version"], "required_version" : "0.9"},\
										"aragorn": { "vpattern" : "ARAGORN v(\d+(\.\d+)+).*", "executable" : ["aragorn"], "version_option" : ["-h"], "required_version" : "1.2.38"},\
										"prodigal": { "vpattern" : "Prodigal V(\d+(\.\d+)+).*", "executable" : ["prodigal"], "version_option" : ["-v"], "required_version" : "2.6.3"}	}
										#todo: paths to excecutables need to be parsed from configs! Here only defaults, assuming executables are in PATH

class version_object(object):
	def __init__(self, version_string):
		templist = [0]*3 #Only comparing major, minor and patch versions. Not accomodating any potential strange special cases here (dependency list is rather small anyway)
		self.version_string = "None"
		if version_string != None:
			if isinstance(version_string, str):
				self.version_string = re.search("\d+(\.\d+)+",version_string).group(0) #this will get rid of prefixes such as "v" but also of non-numeric pre-release designations
				version_digits = version_string.split(".")[:3]
			elif isinstance(version_string, tuple):
				version_digits = version_string[:3]
				self.version_string = ".".join([str(d) for d in version_digits])
			else:
				raise Exception("\nERROR: version_object does not accept input of type '{}'\n".format(type(version_string)))
			for i in range(len(version_digits)):
				templist[i] = int(version_digits[i])
		self.version_tuple = tuple(templist)
		self.major = templist[0]
		self.minor = templist[1]
		self.patch = templist[2]
		
	def __lt__(self, other):
		if isinstance(other, version_object):
			return self.version_tuple < other.version_tuple
		elif isinstance(other, tuple):
			return self.version_tuple < other
		else:
			raise Exception("\nError: version_objects can only be compared against other version_objects or tuples!\n")
	
	def __le__(self, other):
		if isinstance(other, version_object):
			return self.version_tuple <= other.version_tuple
		elif isinstance(other, tuple):
			return self.version_tuple <= other		
		else:
			raise Exception("\nError: version_objects can only be compared against other version_objects or tuples!\n")
				
	def __eq__(self, other):
		if isinstance(other, version_object):
			return self.version_tuple == other.version_tuple
		elif isinstance(other, tuple):
			return self.version_tuple == other	
		else:
			raise Exception("\nError: version_objects can only be compared against other version_objects or tuples!\n")
				
	def __gt__(self, other):
		if isinstance(other, version_object):
			return self.version_tuple > other.version_tuple
		elif isinstance(other, tuple):
			return self.version_tuple > other
		else:
			raise Exception("\nError: version_objects can only be compared against other version_objects or tuples!\n")
			
	def __ge__(self, other):
		if isinstance(other, version_object):
			return self.version_tuple >= other.version_tuple
		elif isinstance(other, tuple):
			return self.version_tuple >= other
		else:
			raise Exception("\nError: version_objects can only be compared against other version_objects or tuples!\n")
			
def get_external_dependency_version_string(toolname):
	cmd = external_dependency_version_dict[toolname]["executable"] + external_dependency_version_dict[toolname]["version_option"]
	versionpattern = external_dependency_version_dict[toolname]["vpattern"]
	try:
		proc = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.STDOUT, text = True, universal_newlines=True)
		loopcounter = 0
		while loopcounter < 10: #expecting version to be mentioned within the first 10 lines...		
			output = proc.stdout.readline().strip()
			versionhit = re.search(versionpattern,output)
			if versionhit:
				return versionhit.group(1)
			loopcounter += 1
	except FileNotFoundError:
		sys.stderr.write("Could not find '{}'!".format(toolname))
		return None
	except OSError:
		sys.stderr.write("Could not execute '{}' for dependency: '{}'!".format(cmd[0], toolname))
		return None	

def get_python_module_version(modulename):
	try:
		m = importlib.import_module(modulename)
	except Exception:
		sys.exit("\nERROR: could not load module '{}'".format(modulename))
	return m.__version__
	
def check_python_dependencies():
	for p in python_dependency_version_dict:
		if python_dependency_version_dict[p]["ist"] < python_dependency_version_dict[p]["soll"]:
			sys.exit("\nERROR: Wrong version of {}! is '{}' but should be '{} or newer!\n".format(p, python_dependency_version_dict[p]["ist"].version_string, python_dependency_version_dict[p]["soll"].version_string))

python_dependency_version_dict = {	"python": {"soll": version_object("3.7"), "ist" : version_object(sys.version_info[:3])},\
									"biopython": {"soll": version_object("1.78"), "ist" : version_object(get_python_module_version("Bio"))} }
		
def check_dependencies(*toolnames, configs):
	check_python_dependencies()
	check_external_dependency(*toolnames, configs=configs)
	
def check_external_dependency(*toolnames, configs):
	sys.stderr.write("\n\tchecking dependencies...\n")
	global external_dependency_version_dict
	for tool in external_dependency_version_dict:
		if tool in configs.settings:
			external_dependency_version_dict[tool]["executable"] = [configs.settings[tool]]
	if len(toolnames) == 0: #by default, check ALL dependencies
		toolnames = [tool for tool in external_dependency_version_dict]
	allisfine = True
	for tool in toolnames:
		sys.stderr.write("\t\t{}...".format(tool))
		isttool = version_object(get_external_dependency_version_string(tool))
		sys.stderr.write(isttool.version_string)
		solltool = version_object(external_dependency_version_dict[tool]["required_version"])
		if isttool < solltool:
			allisfine = False
			sys.stderr.write(" -->WARNING: Incorrect version or excecutable not detected: {}\n".format(tool))
		else:
			sys.stderr.write(" --> OK!\n")
	if not allisfine:
		sys.exit("\n\tERROR: Not all required dependencies for running MDMcleaner are met\n")
	sys.stderr.write("\t-->OK! All dependencies are being met!\n")

# ~ for tool in external_dependency_version_dict:
	# ~ print("{} : {}".format(tool, get_external_dependency_version_string(tool)))
# ~ check_dependency(configs=None)
