'''
This is a general purpose logging utility.  It is intended to record the history
of associated with processing of the data, mainly but not exclusively of
errors that are encountered.
'''

import datetime
import subprocess
import os
import sys

log_file = None

def open_log(file_name, reinitialize=False):
    global log_file
    
    if log_file is not None:
        log_file.close()
    
    if reinitialize:
        mode = 'w'
    else:
        mode = 'a'
    
    log_file = open(file_name, mode)

    commit=get_current_commit()
    if commit is not None:
        log_message(f"Current commit: {commit}")
    else:
        print("Unable to retrieve the current commit.")

def log_message(message):

    if log_file is None:
        sys.exit("Programming Error: Log file is not open. Please call 'open_log()' to open a log file.")


    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_file.write(f'[{timestamp}] {message}\n')
    log_file.flush()
    print(f'[{timestamp}] {message}')

def close_log():
    if log_file is not None:
        log_file.close()


def get_current_commit():

    current_dir=os.getcwd()
    try:
        # Get the absolute path of the current script
        script_path = os.path.abspath(__file__)
        
        # Navigate to the script's directory
        script_directory = os.path.dirname(script_path)
        os.chdir(script_directory)
        
        # Execute 'git rev-parse HEAD' command to get the current commit hash
        output = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
        commit_hash = output.strip().decode('utf-8')
        os.chdir(current_dir)
        
        return commit_hash
    except (subprocess.CalledProcessError, FileNotFoundError):
        os.chdir(current_dir)
        return None


