import os,mmap,sys,subprocess,csv,gzip,pickle,shlex,time
from tempfile import NamedTemporaryFile
from functools import partial
from collections import defaultdict as dd
import numpy as np
import multiprocessing,csv 
cpus = multiprocessing.cpu_count()

mem_bytes = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')  # e.g. 4015976448
mem_mib = mem_bytes/(1024.**2) 
proc_mem = mem_mib / (cpus +1)

allele_dict = dd(str)
for a1,a2 in [('T','A'),('C','G')]:
    allele_dict[a1] = a2
    allele_dict[a2] = a1


def map_alleles(a1,a2):
    """
    Flips alleles to the A strand if neccessary and orders them lexicogaphically
    """
    # check if the A variant is present
    if 'A' not in a1 + a2 and 'a' not in a1+a2:
        # for both/ref and alt map them to the A strand and order each one lexicographically
        a1 = flip_strand(a1)
        a2 = flip_strand(a2)
    # further sorting
    return sorted([a1,a2])


def flip_strand(allele):
    return ''.join([allele_dict[elem.upper()] for elem in allele])

def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise                

def file_exists(fname):
    '''
    Function to pass to type in argparse
    '''
    if os.path.isfile(fname):
        return str(fname)
    else:
        print(fname + ' does not exist')
        sys.exit(1)
        
def mapcount(filename):

    if not os.path.isfile(filename):
        raise ValueError("File doesn't exist")
    
    try:
        return count_lines(filename)
    except:
        return 0
    
def count_lines(filename):
    '''
    Counts line in file
    '''
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    while readline():
        lines += 1
    return lines

def tmp_bash(cmd,check = False):
    

    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")

    os.chmod(scriptFile.name,0o777)
    scriptFile.file.close()

    if check:
        subprocess.check_call(scriptFile.name)
    else:
        subprocess.call(scriptFile.name,stderr = subprocess.DEVNULL)


def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''
   
    file_path,file_root,file_extension = get_path_info(f)

    if 'bgz' in file_extension:
        #print('gzip.open with rb mode')
        open_func = partial(gzip.open, mode = 'rb')
    
    elif 'gz' in file_extension:
        #print('gzip.open with rt mode')
        open_func = partial(gzip.open, mode = 'rt')

    else:
        #print('regular open')
        open_func = open      
    return open_func

def get_path_info(path):
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path,file_root,file_extension

def identify_separator(f):
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(header)
    return dialect.delimiter
    
def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


def basic_iterator(f,separator = None,skiprows = 0,count = False,columns = 'all'):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    '''

    open_func = return_open_func(f)
    if not separator:separator = identify_separator(f)
    i = open_func(f)
    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)          
            yield line
    else:
        row = 0
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            row += 1   
            yield row,line
            
def return_columns(l,columns):
    '''
    Returns all columns, or rather the elements, provided the columns
    '''
    if columns == 'all':
        return l
    elif type(columns) == int:
        return l[columns]
    elif type(columns) == list:
        return list(map(l.__getitem__,columns))

def return_header(f):

    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    delimiter = identify_separator(f)
    header = header.split(delimiter)
    return header
        
        
def fix_header(file_path):
    """
    Strip excess spaces if needed
    """
    #identify separator and fix header if necessary
    header = return_header(file_path)
    header_fix = []
    for elem in header:
        if '  ' in elem:
            while '  ' in elem: elem = elem.replace('  ',' ')
        if ' ' in elem:
            for s in elem.split(' '): header_fix.append(s)
        else:
            header_fix.append(elem)
    return header_fix

       
def isfloat(value):
    try:
        # return True if value is a float and NOT infinite
        value =float(value)
        return ~np.isinf(value)

    except ValueError:
        return False


def pretty_print(string,l = 30):
    l = l-int(len(string)/2)
    print('-'*l + '> ' + string + ' <' + '-'*l)
    

def mapcount_gzip(filename):
    if not os.path.isfile(filename):
        raise ValueError("File doesn't exist")
    try:
        return count_gzip(filename)
    except:
        return 0
    
def count_gzip(myfile):

    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        tmp_bash(f"zcat {myfile} | wc -l > {scriptFile.name} ")
    return int(open(scriptFile.name).read())  

def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, given value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

def load_pos_mapping(chrompos_map):
    '''
    Loads the chrom_pos to ref/alt mapping for finngen
    '''
    out_pickle = chrompos_map + '.pickle'
    if os.path.isfile(out_pickle):
        print('loading chrompos dict -->', out_pickle)

        with open(out_pickle,'rb') as i: pos_dict = pickle.load(i)
    else:
        pos_dict = dd(list)
        print('generating chrompos dict -->', out_pickle)
        iterator = basic_iterator(chrompos_map)
        for entry in iterator:
            chrom_pos,ref,alt = entry
            pos_dict[chrom_pos].append([ref,alt])
        with open(out_pickle,'wb') as o:
            pickle.dump(pos_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
    print('done.')
    return pos_dict

def load_rsid_mapping(rsid_map,inverse = False):
    '''
    Loads the chrompos to rsid mapping
    '''
    out_pickle = os.path.join(rsid_map) + '.pickle'
    if inverse:
        out_pickle += '.chrompos'
    if os.path.isfile(out_pickle):
        print('loading rsid dict -->', out_pickle)
        with open(out_pickle,'rb') as i: rsid_dict = pickle.load(i)
    else:
        print('generating rsid dict --> ',out_pickle)
        rsid_dict = dd(str)
        header = return_header(rsid_map)
        rsid_col = [header.index(elem) for elem in header if 'rs' in elem][0]
        columns = [0,1] if not rsid_col else [1,0]
        if inverse: columns = columns[::-1]
        
        for entry in basic_iterator(rsid_map,columns = columns):
            rsid,chrom_pos = entry
            rsid_dict[rsid] = chrom_pos
        with open(out_pickle,'wb') as o:
            pickle.dump(rsid_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
        print('done.')
    return rsid_dict




def natural_sort(l):
    import re
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def timing_function(some_function):

    """
    Outputs the time a function takes  to execute.
    """

    def wrapper(*args,**kwargs):
        t1 = time.time()
        some_function(*args)
        t2 = time.time()
        print("Time it took to run the function: " + str((t2 - t1)))

    return wrapper



def merge_files(o_file,file_list):
    with open(o_file,'wt') as o:
        for f in file_list:
            with open(f,'rt') as i:
                for line in i:
                    o.write(line)
