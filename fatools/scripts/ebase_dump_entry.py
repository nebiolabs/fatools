#!/usr/bin/python
#

import os
import glob
import zipfile
import fatools.scripts.fa as fa

import psycopg2
import pysftp

def main():

    show_entry = 2782

    try:
        f = open('dbinfo')
        conn_str = f.readlines()[0]
        f.close()
        conn = psycopg2.connect(conn_str)
    except:
        print("I am unable to connect to the database")

    cur = conn.cursor()

    # get all the fields for one entry and print them to the screen
    try:
        ex0 = "SELECT column_name FROM information_schema.columns WHERE table_name = 'ce_experiments' " + \
              "AND table_schema = 'public';"
        cur.execute(ex0)
        fields = [item[0] for item in cur.fetchall()]

        ex1 = "SELECT " + ', '.join(fields) + " FROM ce_experiments WHERE id=" + str(show_entry) + ";"

        cur.execute(ex1)
        rows = cur.fetchall()

        for row in rows:
            for field, value in zip(fields, row):
                print(field, ": ", value)

    except Exception as e:
        print(e)


if __name__ == '__main__':

    main()
