"""
Creates a SMILES batch file with ID's for use with EPI suite and a translation file to translate between a generated id and an InChI
"""
import sqlite3

def db(database = "dataset.db"):
    '''Sets up a database connecion and returns this connection and a cursor object'''
    conn = sqlite3.Connection(database)
    cur = conn.cursor()
    return conn, cur

def add_epi_input(ident, smile):
    '''
    Adds an identifier and smile to the SMILES batch file for EPI suite.
    This file requires to be formatted as [SMILE][space][ID]
    '''

    # Opens or creates file if it doesnt exist yet
    # Will not overwrite any pre existing file!
    with open("epi_input.txt", "a+") as f:
        f.write(f"{smile} {ident}\n")

def add_list(ident, values):
    '''Creates an entry only for each available SMILE in the EPI suite SMILES batch file'''
    
    if values["cas_smile"]:
        add_epi_input(f"{ident}_cs",values["cas_smile"])
        
    if values["cas_cansmile"]:
        add_epi_input(f"{ident}_cc",values["cas_cansmile"])
        
    if values["pc_cansmile"]:
        add_epi_input(f"{ident}_pc",values["pc_cansmile"])
        
    if values["pc_isosmile"]:
        add_epi_input(f"{ident}_pi",values["pc_isosmile"])

def add_translation(ident, inchi):
    '''
    Adds to (and creates if it doesnt exists) a file linking the identification number to a InChI
    Mind that this file is tab separated for readability which is different from the EPI input file!
    '''
    with open("translation.txt", "a+") as f:
        f.write(f"{ident}\t{inchi}\n")

if __name__ == "__main__":
    # sets up database
    conn, cur = db()

    # retrieves InChIs and all available smiles
    cur.execute("""
                SELECT DISTINCT
                    Compound_entries.id,
                    CAS_data.inchi as "CAS inchi",
                    PC_data.inchi as "PC inchi",
                    CAS_data.smile as "CAS smile",
                    CAS_data.canonical_smile as "CAS can_smile",
                    PC_data.canonical_smiles as "PC_can_smile",
                    PC_data.isometric_smiles as "PC_iso_smile"
                FROM
                    Compound_entries
                LEFT JOIN
                    CAS_data ON CAS_data.cas_rn = Compound_entries.CAS_data_id
                LEFT JOIN
                    PC_data ON PC_data.cid = Compound_entries.PC_data_id
                WHERE
                    NOT (Compound_entries.CAS_data_id IS NULL AND Compound_entries.PC_data_id IS NULL)
                GROUP BY
                    Compound_entries.id
                """)
    results = cur.fetchall()

    # Sets up a dictionary to store InChI as keys and smiles as data
    inchi_list = {}

    # Process results
    for i in results:
        # Takes InChI to use as the InChI provided by the CAS register if available, otherwise uses the InChI to the one from PubChem
        # Due to the design of the identification script these should be the same if both are available
        cas_inchi = i[1]
        pc_inchi = i[2]
        inchi = cas_inchi if cas_inchi else pc_inchi

        # Only adds InChI unique InChI to the list to process
        if not inchi in inchi_list.keys():
            inchi_list[inchi] = {"cas_smile":i[3],"cas_cansmile":i[4],"pc_cansmile":i[5],"pc_isosmile":i[6]}

    # Create both output files
    for index, i in enumerate(inchi_list.keys()):
        add_translation(index, i)
        add_list(index, inchi_list[i])

    # Prints some statistics to indicate the script has finished 
    print(f"Created output files for {len(inchi_list)} compounds")
