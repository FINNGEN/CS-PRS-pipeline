#!/usr/bin/env python3
import shlex,os
from subprocess import  call,check_output
import argparse,datetime

curr_path =os.path.dirname(os.path.realpath(__file__))
root_path = '/'.join(curr_path.split('/')[:-1])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build Docker file for variant filtering")

    parser.add_argument("--docker", type= str,
                        help="path of docker file",default = curr_path + '/Dockerfile')
    
    parser.add_argument("--image", type= str,
                        help="name of image",default = 'cs-prs')
    parser.add_argument("--version", type= str,
                        help="version value, e.g.0.001",required = True)
    parser.add_argument("--push",action = 'store_true')
    parser.add_argument("--args",type = str,default = '')
    args = parser.parse_args()
    
    docker_path = "eu.gcr.io/finngen-refinery-dev/"
    
    cmd = f"docker build -t {docker_path}{args.image}:{args.version} -f {args.docker} {root_path} {args.args}"
    print(cmd)
    call(shlex.split(cmd))

    if args.push:
        current_date = datetime.datetime.today().strftime('%Y-%m-%d')
        cmd = f"gcloud docker -- push {docker_path}{args.image}:{args.version}"
        with open('./docker.log','a') as o:o.write(' '.join([current_date,cmd]) + '\n')
        call(shlex.split(cmd))
