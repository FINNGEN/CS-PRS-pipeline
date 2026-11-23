#!/usr/bin/python3
import shlex, os, argparse, datetime, subprocess
from subprocess import Popen, PIPE, call, check_output
from pathlib import Path

path = Path(__file__).parent.absolute()

# Path Constants
REFINERY_PATH = "europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/"
SB_PATH = "eu.gcr.io/finngen-sandbox-v3-containers/"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build Docker file for CS-PRS-pipeline")

    parser.add_argument("--image", type=str, help="name of image", default='cs-prs')
    parser.add_argument("--version", type=str, help="version value, e.g. 0.001", required=True)
    parser.add_argument("--push", action='store_true')
    parser.add_argument("--args", type=str, default='')

    # New argument for choosing the registry path
    parser.add_argument("--registry",
                        choices=['refinery', 'sandbox'],
                        default='refinery',
                        help="Choose the destination registry (default: refinery)")

    args = parser.parse_args()

    # Determine which base path to use
    base_path = REFINERY_PATH if args.registry == 'refinery' else SB_PATH
    full_image_tag = f"{base_path}{args.image}:{args.version}"

    # Build Command
    build_cmd = f"docker build -t {full_image_tag} -f {os.path.join(path, 'Dockerfile')} {path.parent} {args.args}"
    print(f"Running: {build_cmd}")
    call(shlex.split(build_cmd))

    # Push Command
    if args.push:
        push_cmd = f"docker push {full_image_tag}"
        print(f"Running: {push_cmd}")
        call(shlex.split(push_cmd))
