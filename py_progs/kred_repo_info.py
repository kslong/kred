#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Get information about the date and commit in the 
git repository of the sofware being run


Command line usage (if any):

    usage: get_kred_repo_infoa


Description:  

Primary routines:

    get_kred_repo_info

Notes:

    This is just a simple routine so that any time we
    can obtain the data and time of the version of
    kred sofware used to run a file
                                       
History:

240501 ksl Coding begun

'''



    
import subprocess
import os



def get_kred_git_commit_info():
    repo_path=git_repo_path = os.environ.get('KRED')
    try:
        # Run git command to get the latest commit hash and commit date
        result_hash = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=repo_path, stdout=subprocess.PIPE)
        commit_hash = result_hash.stdout.strip().decode('utf-8')

        result_date = subprocess.run(['git', 'show', '-s', '--format=%ci', commit_hash], cwd=repo_path, stdout=subprocess.PIPE)
        commit_date = result_date.stdout.strip().decode('utf-8')

        return commit_hash, commit_date
    except Exception as e:
        print(f"Error getting Git commit information: {e}")
        return None, None


# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
    import sys
    if len(sys.argv)==1:
        print(get_kred_git_commit_info())
    else:
        print (__doc__)
