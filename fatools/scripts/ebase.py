#!/usr/bin/python
#

import os
import glob
import zipfile
import fatools.scripts.fa as fa
import time

import psycopg2
import pysftp

class MyPeak():

    def __init__(self, dye, size_s, size_bp, area_s, area_bp, height):

        self.dye     = dye
        self.size_s  = size_s
        self.size_bp = size_bp
        self.area_s  = area_s
        self.area_bp = area_bp
        self.height  = height

    def __repr__(self):
        return "<P: dye %3s | size_s %3d | size_bp %4d | area_s %3d | area_bp %4d | height %5d >" % (
            self.dye, self.size_s, self.size_bp, self.area_s, self.area_bp, self.height )

def main():

    start_time = time.time()
    
    nrecords = 1
    use_db = True
    #trace_dir = "1_and_2_base_extensions" 
    #trace_dir = "single_base_extension" 
    #trace_dir = "trace-Mar-22-2010-Aug13-14-23-50-mismatch5" #"trace-Aug14-2012_primase_100-400nM_CE13531"
    trace_dir = ""
    file_list = ""#"Blank.fsa"
    scp_files = True # only used if use_db = True
    overwrite_output = True
    
    if scp_files:
        f = open('dbinfo')
        lines = f.readlines()
        f.close()
        username = lines[1].rstrip('\n')
        keyfile = lines[2].rstrip('\n')

        sftp = pysftp.Connection('ebase-c.neb.com', username=username, private_key=keyfile)

    files = {} # dictionary containing directories and names of files within directories

    if use_db:

        conn = None

        # open the database
        try:
            f = open('dbinfo')
            conn_str = f.readlines()[0]
            f.close()
            conn = psycopg2.connect(conn_str)
        except:
            print("I am unable to connect to the database")

        cur = conn.cursor()

        # get sets of filenames and directories
        """ex = "SELECT id, submission_date, data_file_name, peaks_table_file_name " + \
             "FROM ce_experiments " + \
             "WHERE submission_date > '1/1/2016' " + \
             #"ORDER BY submission_date DESC " + \
             "LIMIT " + str(nrecords) + " " + \
             "OFFSET 20 ";"""

        ex = "SELECT id, submission_date, data_file_name, peaks_table_file_name " + \
             "FROM ce_experiments " + \
             "LIMIT " + str(nrecords) + " " + \
             "OFFSET 20 ";

        cur.execute(ex)
        rows = cur.fetchall()

        # go to output directory
        if not os.path.isdir("../output"):
            os.makedirs("../output")

        done = False
        for row in rows:

            if done: break

            id, submission_date, zipped_files, peaks_table_file_name = row[0], row[1], row[2], row[3]
            file_root = zipped_files[:-4]

            print("submission date: ", submission_date)

            if file_root == "Archive": # skip these for now... the script doesn't like multiple directories with the same name
                continue

            if trace_dir!="":
                if file_root == trace_dir:
                    done=True
                else:
                    continue

            basedir = "/var/www/ebase/shared/shared_uploads/ce_experiments/" + str(id)
            print("analyzing ", basedir+"/"+file_root)

            full_zipped_files_name = basedir + "/" + zipped_files
            full_peaks_table_file_name = basedir + "/" + peaks_table_file_name

            if scp_files:
                sftp.get(full_zipped_files_name)
                full_zipped_files_name = zipped_files
                sftp.get(full_peaks_table_file_name)
                    

            if not os.path.isdir(file_root):
                os.makedirs(file_root)
                os.chdir(file_root)
                os.rename("../"+zipped_files,zipped_files)
                os.rename("../"+peaks_table_file_name,peaks_table_file_name)
                    
                with zipfile.ZipFile(full_zipped_files_name) as zip_ref:
                    zip_ref.extractall(".")
                if scp_files:
                    os.remove(zipped_files)
                    
                # sometimes there is another directory with the same name
                # containing all the fsa files
                if os.path.isdir(file_root):
                   os.chdir(file_root)
                   fsa_files = glob.glob("*fsa")
                   for f in fsa_files:
                       os.rename(f,"../"+f)
                   os.chdir("..")
                os.chdir("..")
            files[id] = file_root

        if scp_files:
            sftp.close()
                
    # read from local files instead of database
    else:
        files['2'] = trace_dir

    for key, root_dir in files.items():

        fa_args = ['--align',
                   '--panel=GS120LIZ',
                   '--ladder_rfu_threshold=0.25',
                   '--nonladder_rfu_threshold=150',
                   '--nonladder_peak_window=5',
                   '--call',
                   #'--ladderplot',
                   '--allelemethod=localsouthern',
                   '--baselinemethod=minimum',
                   '--baselinewindow=51',
                   '--listpeaks',
                   '--peaks_format=peakscanner',
                   '--verbose=0',
                   #'--plot=B,O',
                   #'--range=41,47',
                   #'--range=-15,135',
                   '--indir='+root_dir]

        if file_list!="":
            fa_args.append('--file='+file_list)

        if file_list!="" and len(file_list.split(','))==1:
            fa_args.append('--outfile='+root_dir+"/"+root_dir+"_"+file_list+".out")
        else:
            fa_args.append('--outfile='+root_dir+"/"+root_dir+".out")

        print("fa_args: ", fa_args)

        parser = fa.init_argparser()
        args = parser.parse_args(fa_args)
        fa.main(args)

        # do some cleanup of the output directory
        if use_db:

            os.chdir(root_dir)

            # first get list of badfiles
            f = open(root_dir+'_badfiles.out','r')
            badfiles = f.read().splitlines()
            f.close()

            badfiles = [ i.split(':')[1].strip() for i in badfiles ]
            
            for fsa_filename in glob.glob("*.fsa"):
                if fsa_filename not in badfiles:
                    os.remove(fsa_filename)
                else:
                    print("keeping bad file ",fsa_filename)

            # copy ebase.py and params.py to output directory
            print("cwd: ", os.getcwd())
            import shutil
            shutil.copy("../ebase.py",".")
            shutil.copy("../../lib/params.py", ".")
            
            os.chdir("..")

            if os.path.isdir("../output/"+root_dir):
                if overwrite_output: # overwrite output directory
                    import shutil
                    shutil.rmtree("../output/"+root_dir)
                    os.rename(root_dir,"../output/"+root_dir)
                else: # append a number to the output directory
                    output_dirs = glob.glob("../output/"+root_dir+"_*")
                    n_existing_dirs = 1+len(output_dirs)
                    os.rename(root_dir,"../output/"+root_dir+"_"+
                              str(n_existing_dirs))
                    
            else:
                os.rename(root_dir,"../output/"+root_dir)

    elapsed_time = time.time() - start_time
    print("elapsed time: ", elapsed_time)
    
if __name__ == '__main__':

    main()
