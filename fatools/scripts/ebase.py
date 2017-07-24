#/usr/bin/python
#

import os, glob, zipfile
import psycopg2
import argparse
import importlib
import fa
import pysftp

class MyPeak ():

    def __init__(self, dye, size_s, size_bp, height):

        self.dye     = dye
        self.size_s  = size_s
        self.size_bp = size_bp
        self.height  = height

    def __repr__(self):
        return "<P: dye %3s | size_s %3d | size_bp %4d | height %5d >" % (
            self.dye, self.size_s, self.size_bp, self.height )

def main():
    
    nrecords = 2
    print_all = False
    use_db = True
    scp_files = True # only used if use_db = True

    if scp_files:
        f = open('dbinfo')
        lines = f.readlines()
        f.close()
        username = lines[1].rstrip('\n')
        keyfile = lines[2].rstrip('\n')

        sftp = pysftp.Connection('ebase-c.neb.com', username=username, private_key=keyfile)

    files = {} # dictionary containing directories and names of files within directories

    if use_db:

        # open the database
        try:
            f = open('dbinfo')
            conn_str = f.readlines()[0]
            f.close()
            conn = psycopg2.connect(conn_str)
        except:
            print("I am unable to connect to the database")

        cur = conn.cursor()

        # get all the fields for one entry and print them to the screen
        if print_all:
            try:
                ex0 = "SELECT column_name FROM information_schema.columns WHERE table_name = 'ce_experiments' AND table_schema = 'public';"
                cur.execute(ex0)
                fields = [item[0] for item in cur.fetchall()]

                ex1 = "SELECT " + ', '.join(fields) + " FROM ce_experiments LIMIT " + str(nrecords) + ";"
                cur.execute(ex1)
                rows = cur.fetchall()

                for row in rows:
                    for field, value in zip(fields, row):
                        print(field, ": ", value)

            except Exception as e:
                print(e)

        # get sets of filenames and directories
        try:
            ex = "SELECT id, data_file_name, peaks_table_file_name FROM ce_experiments LIMIT " + str(nrecords) + ";"
            cur.execute(ex)
            rows = cur.fetchall()

            for row in rows:
                
                id, zipped_files, peaks_table_file_name = row[0], row[1], row[2]

                print("peaks_table_file_name: ",peaks_table_file_name)

                basedir = "/var/www/ebase/shared/shared_uploads/ce_experiments/" + str(id)
                full_zipped_files_name = basedir + "/" + zipped_files
                full_peaks_table_file_name = basedir + "/" + peaks_table_file_name

                if scp_files:
                    sftp.get(full_zipped_files_name)
                    full_zipped_files_name = zipped_files
                    sftp.get(full_peaks_table_file_name)
                    
                fileroot = zipped_files[:-4]
                if not os.path.isdir(fileroot):
                    os.makedirs(fileroot)
                    os.chdir(fileroot)
                    os.rename("../"+zipped_files,zipped_files)
                    os.rename("../"+peaks_table_file_name,peaks_table_file_name)
                    
                    with zipfile.ZipFile(full_zipped_files_name) as zip_ref:
                        zip_ref.extractall(".")
                    if scp_files:
                        os.remove(zipped_files)
                    os.chdir("..")
        
                file_list = ""
                files[id] = (fileroot, peaks_table_file_name, file_list) 

            if scp_files:
                sftp.close()
                
        except Exception as e:
            print(e)

    # read from local files instead of database
    else:
        fileroot = "trace-Jul-25-2013-Jul26-13-32-28_full"
        file_list = "GL8-044D3-10.fsa"
        peaks_table_file_name = "GL8-044.csv"

        files['2'] = (fileroot, peaks_table_file_name, file_list)


    for key, info in files.items():
        
        file_list = info[2]
        file_root = info[0]

        fa_args = ['--align',
                   '--panel=GS120LIZ',
                   '--ladder_rfu_threshold=1000',
                   '--nonladder_rfu_threshold=50',
                   '--call',
                   #'--plot', '--ladderplot',
                   '--allelemethod=leastsquare',
                   '--baselinemethod=median',
                   '--baselinewindow=399',
                   '--listpeaks',
                   '--peaks_format=peakscanner',
                   '--verbose=0',
                   '--indir='+file_root]

        if file_list!="":
            fa_args.append('--file='+file_list)

        if file_list!="" and len(file_list.split(','))==1:
            fa_args.append('--outfile='+fileroot+"/"+fileroot+"_"+file_list+".out")
        else:
            fa_args.append('--outfile='+fileroot+"/"+fileroot+".out")

        try:
            parser = fa.init_argparser()
            args = parser.parse_args(fa_args)
            fa.main(args)

        except Exception as e:
            print(e)

        # do some cleanup of the output directory
        if use_db:

            os.chdir(file_root)

            # first get list of badfiles
            f = open(file_root+'_badfiles.out','r')
            badfiles = f.read().splitlines()
            f.close()

            badfiles = [ i.split(':')[1].strip() for i in badfiles ]
            
            print("badfiles: ", badfiles)
            
            for fsa_filename in glob.glob("*.fsa"):
                #print("file: ", fsa_filename)
                if fsa_filename not in badfiles:
                    os.remove(fsa_filename)
                else:
                    print("keeping file ",fsa_filename)
            os.chdir("..")
        

if __name__ == '__main__':

    main()
