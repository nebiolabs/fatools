#/usr/bin/python
#

import os, glob, zipfile
import psycopg2
import argparse
import importlib
import fa

# try to connect

nrecords = 1
do_all = False

#fields = ['id','name','submission_date','number_of_wells','operator_id','data_file_name',
#          'data_content_type','data_file_size','data_updated_at','peaks_table_file_name',
#          'peaks_table_content_type','peaks_table_updated_at','peaks_table_file_size',
#          'experiment_type','well_set_id','notes','state','excluded','exclusion_reason',
#          'created_at','updated_at','creator_id','peak_th','peak_num','peak_dye',
#          'ce_bin_set_id','bin_selection','flag_th']
          
try:
    f = open('dbinfo')
    conn_str = f.readlines()[0]
    conn = psycopg2.connect(conn_str)
except:
    print("I am unable to connect to the database")

cur = conn.cursor()

# import fa 
#try:
#    M = importlib.import_module('fa')

#except ImportError:
#    print('Cannot import script name: fatools.scripts.fa')
#    raise

try:
    if do_all:
        for field in fields: 

            if field != 'data_file_name': continue
        
            ex = "SELECT "+field+" FROM ce_experiments LIMIT " + str(nrecords) + ";"
            cur.execute(ex)
            rows = cur.fetchall()

            for row in rows:
                for entry in row:
                    print(field,": ",entry)

    else:
        ex = "SELECT id, data_file_name, peaks_table_file_name FROM ce_experiments LIMIT " + str(nrecords) + ";"
        cur.execute(ex)
        rows = cur.fetchall()
        for row in rows:

            id, zipped_files, peaks_files = row[0], row[1], row[2]
    
            print("peaks_files: ",peaks_files)

            full_zipped_files_name = "/var/www/ebase/shared/shared_uploads/ce_experiments/" + str(id) + "/" + zipped_files

            fileroot = zipped_files[:-4]
            if not os.path.isdir(fileroot):
                os.makedirs(fileroot)
                os.chdir(fileroot)
                with zipfile.ZipFile(full_zipped_files_name) as zip_ref:
                    zip_ref.extractall(".")

            os.chdir(fileroot)
            file_list = ""
            for file in glob.glob("*.fsa"):
                if file_list == "":
                    file_list = fileroot+"/"+file
                else:
                    file_list = file_list + ", " + fileroot+"/"+file
            print("file_list: ",file_list)
            os.chdir("..")

            fa_args = ['--align',
                       '--panel=GS120LIZ',
                       '--ladder_rfu_threshold=1000',
                       '--nonladder_rfu_threshold=100',
                       '--call',
                       '--allelemethod=leastsquare',
                       '--listpeaks',
                       '--verbose=7',
                       '--indir='+fileroot,
                       '--file='+file_list,
                       '--outfile='+fileroot+"/"+fileroot+'.out']
            
            parser = fa.init_argparser()
            args = parser.parse_args(fa_args)
            fa.main(args)

except Exception as e:
    print(e)


