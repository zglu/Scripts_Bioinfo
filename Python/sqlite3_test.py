#! /usr/bin/env python
# -*- coding: utf-8 -*-  

import sqlite3

# establish a connection
connection = sqlite3.connect("/Documents/PersonDb.db")
# hold a cursor
cursor = connection.cursor()

# add new data
new_person = [
    ("Double-decker", "Bus", 2008),
    ("Electric", "Bike", 2006),
]

for p in new_person:
    format_str = """INSERT INTO Person (firstName, lastName, birthYear)
        VALUES ("{first}", "{last}", "{birthyear}");"""
    sql_command = format_str.format(first=p[0], last=p[1], birthyear=p[2])
    cursor.execute(sql_command)

## remove duplicates
cursor.execute("DELETE FROM Person WHERE rowid NOT IN (SELECT min(rowid) FROM Person GROUP BY firstName, lastName, birthYear);")

## Fetch all remaining records
cursor.execute("SELECT * FROM Person")
print("Fetchall:")
all_records = cursor.fetchall()
for r in all_records:
    print(r)

## Fetch the next record
cursor.execute("SELECT * FROM Person")
print("\nFetch one:")
one_record = cursor.fetchone()
print(one_record)

## close cursor
cursor.close()

# Discard changes
# connection.rollback()

# Save the changes and close the connection:
connection.commit()
connection.close()
