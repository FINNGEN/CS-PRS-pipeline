import os,mmap,sys,subprocess,csv,gzip,pickle
from tempfile import NamedTemporaryFile
from functools import partial
from collections import defaultdict as dd

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
    float(value)
    return True
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

    i = 0
    with gzip.open(myfile, 'rb') as f:
        for i, l in enumerate(f,1):
            pass
    return i


def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, given value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

def load_pos_mapping(rsid_map,out_path):
    '''
    Loads the chrom_pos to ref/alt mapping for finngen
    '''
    out_pickle = os.path.join(out_path,'chrompos.pickle')
    if os.path.isfile(out_pickle):
        with open(out_pickle,'rb') as i: pos_dict = pickle.load(i)
    else:
        pos_dict = dd(list)
        print('generating pickle dict...')
        iterator = basic_iterator(chrompos_map)
        for entry in iterator:
            chrom_pos,ref,alt = entry
            pos_dict[chrom_pos].append([ref,alt])
        with open(out_pickle,'wb') as o:
            pickle.dump(pos_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
           
    return pos_dict

def load_rsid_mapping(rsid_map,out_path):
    '''
    Loads the chrompos to rsid mapping
    '''
    out_pickle  = os.path.join(out_path,'rsid.pickle')     
    if os.path.isfile(out_pickle):
        with open(out_pickle,'rb') as i: rsid_dict = pickle.load(i)
    else:
        print('generating rsid dict...')
        rsid_dict = dd(str)
        header = return_header(rsid_map)
        rsid_col = [header.index(elem) for elem in header if 'rs' in elem][0]
        columns = [0,1] if not rsid_col else [1,0]
        for entry in basic_iterator(rsid_map,columns = columns):
            rsid,chrom_pos = entry
            rsid_dict[rsid] = chrom_pos
        with open(out_pickle,'wb') as o:
            pickle.dump(rsid_dict,o,protocol = pickle.HIGHEST_PROTOCOL)
        print('done.')
    return rsid_dict
