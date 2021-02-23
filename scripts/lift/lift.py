#! /usr/bin/env python3

import argparse,gzip,sys,subprocess,os,shlex
from functools import partial
from tempfile import NamedTemporaryFile


def get_dat_var(line, index,sep):
    d = line[index].split(sep)
    if len(d)<4:
        print("WARNING: Not properly formatted variant id in line: " + line, file=sys.stderr )
        return None
    return d


def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''

    file_path = os.path.dirname(f)
    basename = os.path.basename(f)
    file_root, file_extension = os.path.splitext(basename)

    result = subprocess.run(['file', f], stdout=subprocess.PIPE)
    gzip_bool = True if 'gzip' in result.stdout.decode("utf-8") else False

    if gzip_bool:
        open_func = partial(gzip.open, mode = 'rt')

    else:
        open_func = open
        
    return open_func

def lift(args):
    '''
    Finds the proper header column names and calls the liftover.sh script
    '''
    open_func = return_open_func(args.file)
    with open_func(args.file) as res:
        #skip first line anyways
        header = res.readline().rstrip("\n").split(args.sep)
        print(args.info)
        print(header)
        if args.var:
            if not args.numerical:
                args.var = header.index(args.var)

            joinsortargs =f"--var {args.var+1}"
            get_dat_func = partial(get_dat_var,index=args.var,sep = args.var_sep)

        elif args.info:
            if not args.numerical:
                if not all(elem in header for elem in args.info):
                    # check that all columns are present
                    raise Exception(f"Given columns not in data. Missing: { [ elem for elem in args.info if elem not in header] }")

                args.info = [header.index(elem) for elem in args.info ]
                
            chrom,pos,ref,alt = args.info
            joinsortargs = f"--chr {chrom+1} --pos {pos+1} --ref {ref+1} --alt {alt+1}"
            get_dat_func = lambda line:  (line[chrom], line[pos], line[ref], line[alt])

        tmp_bed = NamedTemporaryFile(delete=True)
        with open(tmp_bed.name, 'w') as out:
            for line in res:
                vardat = get_dat_func(line.strip().split())
                string = "{}\t{}\t{}\t{}".format("chr"+vardat[0], str(int(vardat[1])-1), str(int(vardat[1]) + max(len(vardat[2]),len(vardat[3])) -1), ":".join([vardat[0],vardat[1],vardat[2],vardat[3]])) + "\n"
                out.write(string)

    #change working dir to args.out so i don't have to move errors and variants_lifted
    os.chdir(args.out)
    cmd = f"{args.scripts_path}liftOver  {tmp_bed.name} {args.chainfile} variants_lifted errors"
    print(cmd)
    subprocess.run(shlex.split(cmd))
    with open('errors', 'r') as errs:
        err_count=0
        for l in errs:
            err_count+=1

        if(err_count>0):
            print(f'WARNING. {int(err_count/2)} variants could not be correctly lifter. Consult file {args.out+"/" if args.out else ""}errors for details')

    #if args.sep:
        # easyoptions seem not to be able handle 5 switched ?!?!? not implemented
        #args.sep=args.sep.replace(" ","space")
        #joinsortargs=f'--sep {args.sep} {joinsortargs}'

    joinsort = f"{os.path.join(args.scripts_path,'joinsort.sh')}"
    ##subprocess.run(shlex.split(f"chmod +x {joinsort}"))
    joincmd = f"{joinsort} {args.file} variants_lifted {joinsortargs}"
    subprocess.run(shlex.split(joincmd))
   

if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Add liftover positions to summary file')
    parser.add_argument("file", help=" Whitespace separated file with either single column giving variant ID in chr:pos:ref:alt or those columns separately")
    parser.add_argument("--chainfile", help=" Chain file for liftover",required = True)
    parser.add_argument('-o',"--out", help=" Folder where to save the results",default = os.getcwd())
    parser.add_argument("--sep", help="column separator in file to be lifted. Default tab", default='\t', required=False)
    group = parser.add_mutually_exclusive_group(required = True)
    group.add_argument('--info',nargs =4, metavar = ('chr','pos','ref','alt'), help = 'Name of columns')
    group.add_argument("--var",nargs = 2,help ="Column name of snpid and separator",metavar = ('snpid','snp_sep'))

    args = parser.parse_args()

    args.file = os.path.abspath(args.file)
    args.chainfile = os.path.abspath(args.chainfile)
    args.out = os.path.abspath(args.out)
    
    if args.var:
        args.var,args.var_sep = args.var
    # checks if var/info are numerical or strings
    args.numerical = False
    if args.var and args.var.isdigit():
        args.var = int(args.var)
        args.numerical = True
    if args.info and  all(elem.isdigit() for elem in args.info):
        args.info = list(map(int,args.info))
        args.numerical = True

    args.scripts_path = '/'.join(os.path.realpath(__file__).split('/')[:-1]) + '/'
    args.file = os.path.abspath(args.file)
    lift(args)
