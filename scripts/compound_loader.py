"""Inserts compounds and related experiment ID from a text file into a database file."""
import sqlite3

if __name__ == "__main__":
    # Opens textfile containing [compound_name][tab][experiment_id] formatted data
    with open("compound_entries.txt") as file:
        entries = file.read().splitlines()

    # Creates a database connection
    conn = sqlite3.connect("dataset.db")
    cur = conn.cursor()

    # Sets up a list of data to insert
    to_insert = []

    for entry in entries:
        # Split line on tab giving a compound name at index 0 and an experiment id at index 1
        entry = entry.split("\t")

        # Tuples of data are prepared by removing any whitespace characters that might be artifacts from the importing process
        to_insert.append((entry[0].strip(),entry[1].strip()))

    # Insert all into the the database, then commit the changes and close the connection
    cur.executemany("INSERT INTO Compound_entries(compound_name, experiment_id) VALUES (?, ?)",to_insert)
    conn.commit()
    conn.close()
