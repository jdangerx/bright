#!/usr/bin/env python

############################
#### Standard Libraries ####
############################
from __future__ import print_function
import os
import sys
import time
import shutil
import logging
import subprocess

from optparse import OptionParser

###########################
### Extension Libraries ###
###########################
import tables as tb

##########################
#### Custom Libraries ####
##########################
import isoname
import metasci
import metasci.nuke as msn
import metasci.graph as msg

from metasci.colortext import failure

#################
### CHAR Libs ###
#################
import envchar

from run.pbs import Pbs
from run.bash import Bash

run_switch = {'': Bash, 
              'BASH': Bash,
              'bash': Bash,
              'PBS': Pbs,
              'pbs': Pbs,
              'Torque': Pbs,
              'torque': Pbs,
              }

from n_code_serpent import NCodeSerpent

n_code_switch = {'': NCodeSerpent, 
                 'sss': NCodeSerpent, 
                 'Serpent': NCodeSerpent, 
                 'serpent': NCodeSerpent, 
                 }

def main():
    global defchar

    ###########################
    ### Command Line Parser ###
    ###########################
    usage = "usage: %prog [options] confchar"
    parser = OptionParser(usage)

    parser.add_option("-v", "--verbose", action="store_true", dest="VERBOSE", 
        default=False, help="Gives extra info while running.")

    parser.add_option("-i", "--input", action="store_true", dest="MAKE_INPUT", 
        help="Makes the transport calculation input deck.")

    parser.add_option("-r", "--run", action="store_true", dest="RUN_TRANSPORT", 
        default=False, help="Run the transport calculation.")

    parser.add_option("-d", "--dry-run", action="store_false", dest="RUN_TRANSPORT", 
        help="Dry Run. Do NOT run the transport calculation.")

    parser.add_option("-a", "--analyze", action="store_true", dest="RUN_ANALYSIS", 
        default=False, help="Run analysis on database.")

    parser.add_option("-b", "--burnup",  action="store_true", dest="RUN_BURNUP", 
        default=False, help="Run the burnup calculation.")

    parser.add_option("-x", "--xs",  action="store_true", dest="RUN_XS_GEN", 
        default=False, help="Run the cross-section generation calculation.")

    parser.add_option("-m", "--delta-mass",  action="store_true", dest="RUN_DELTAM", 
        default=False, help="Run the initial isotope sensitivity calculation.")

    parser.add_option("-c", "--clean", action="store_true", dest="CLEAN", 
        help="Cleans the reactor direactory of current files.")

    parser.add_option("-l", "--local", action="store_true", dest="LOCAL", 
        default=True, help="Run or Fetch files locally.")

    parser.add_option("-s", "--server", action="store_false", dest="LOCAL", 
        help="Run or Fetch files from a remote server.")

    parser.add_option("-f", "--fetch", action="store_true", dest="FETCH_FILES", default=False, 
        help="Fetches files from the remote server. Does not run transport, even if -r is set. Automatically sets -s.")

    parser.add_option("-p", "--pid", action="store_true", dest="PID", default=False, 
        help="Finds the process identification number of a current transport run. Sets -d.")

    parser.add_option("-k", "--kill", action="store_true", dest="KILL_TRANSPORT", 
        default=False, help="Kills the current transport run. Sets -p.")

    parser.add_option("--cwd", action="store_true", dest="CWD", default=False, 
        help="Run char in the current working directory.")

    parser.add_option("--ui", action="store_true", dest="UI", default=False, 
        help="Launches the char ui.")

    (options, args) = parser.parse_args()

    # Try launching ui before anything else
    if options.UI:
        # Test to see if ui library is installed
        try:
            from .ui import app
        except ImportError:
            print(failure("Please install the Enthought Tool Suite (ETS) for CHAR UI."))
            raise SystemExit

        # Open UI
        application = app.Application()
        #application.rx_h5_path = "/home/scopatz/MultiGroupPaper/DataXS/lwr/lwr.h5"
        application.configure_traits()

        # Clean-up UI
        if application.rx_h5 is not None:
            application.rx_h5.close()

        raise SystemExit

    # Make sure we have a configureation file before proceeding
    if len(args) == 0:
        print(failure("Please specify a configuration file for CHAR."))
        raise SystemExit

    # Load the CHAR definition file into its own env namespace
    env = {}
    absolute_path = os.path.abspath(args[0])
    execfile(absolute_path, {}, env)

    # Add command line arguments to env
    env['options'] = options
    env['args'] = args

    # Update defchar adding more useful values.
    env = envchar.update_env(env)

    #intial command-line options protocol.
    if options.KILL_TRANSPORT:
        options.PID = True                      #Ensures that the PID is found in order that it mak be killed.

    if options.PID:
        options.RUN_XS_GEN = False
        options.RUN_DELTAM = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.

    if options.FETCH_FILES:
        options.RUN_XS_GEN = False
        options.RUN_DELTAM = False
        options.RUN_TRANSPORT = False            #Ensures that transport calculation is not initiated while fetching files.
        options.LOCAL = False                   #Ensures that ssh package is loaded.

    ################
    #### Script ####
    ################

    # Prep work
    if not options.CWD:
        if env['options'].CLEAN:
            metasci.safe_remove(env['reactor'], True)

        if env['reactor'] not in os.listdir('.'):
            os.mkdir(env['reactor'])

        os.chdir(env['reactor'])
        shutil.copyfile(absolute_path, 'defchar.py')

    # Start up logger
    logger = logging.getLogger('char')
    hdlr = logging.FileHandler('char.log')
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.INFO)
    env['logger'] = logger

    # Set the transport code
    n_coder = n_code_switch[env['transporter']]
    env['n_code'] = n_coder(env)

    # Get the run controller
    runner = run_switch[env['scheduler']]
    runchar = runner(n_code, env)

    # Make the input file unless otherwise specified.
    if (options.MAKE_INPUT) and (not options.FETCH_FILES) and (not options.PID):
        n_code.make_input()

    # Check a bunch of run conditions
    if options.RUN_TRANSPORT:
        # Run Transport code
        runchar.make_run_script(n_code)

        if options.LOCAL:
            runchar.run_locally()
        else:
            runchar.run_remotely()

    elif options.RUN_ANALYSIS:
        n_code.analyze_deltam()

    elif options.RUN_BURNUP or options.RUN_XS_GEN or options.RUN_DELTAM:
        # Make tranumatrion libraries by executing the as a separate step from 
        # the cross-section generation
        if options.RUN_BURNUP:
            n_code.run_burnup()

        # Make Cross-sections as a separate step from the burnup calculation
        if options.RUN_XS_GEN:
            n_code.run_xs_gen()

        # Run initial isotope sensitivity calculation
        if options.RUN_DELTAM:
            n_code.run_deltam()

    elif options.FETCH_FILES:
        # Fetches files from remote server
        runchar.fetch()
    elif options.PID:
        # Finds and prints the PID of CHAR
        runchar.pid()
    elif options.KILL_TRANSPORT:
        # Finds and kills CHAR
        runchar.kill()

    # Clean up
    if not options.CWD:
        os.chdir('..')


# Run CHAR
if __name__ == '__main__':
    main()
