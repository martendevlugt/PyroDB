import sqlite3
import os
import cas_api
import pubchempy
import json

def db(db_path:str = "dataset.db"):
    '''
    Establish DB connection
    '''
    
    conn = sqlite3.Connection(db_path)
    cur = conn.cursor()
    
    return conn, cur

def add_pc_data(pc_data, conn, cur):
    '''
    Add pubchem data if this cid is not yet in the dataset
    '''

    # Creates dict version from the PubChemPy Compound object
    pc_data = pc_data.to_dict()
    
    # Check for existing records for this identifier in the database, a new record is only added if it does not exist yet.
    cur.execute('SELECT * FROM PC_data WHERE cid = ?', (pc_data['cid'],) )
    results = cur.fetchall()

    if not results:
        # Create a tuple of the data to store in the correct order
        write_data = (pc_data['cid'],
                      str(pc_data['elements']),
                      str(pc_data['atoms']),
                      str(pc_data['bonds']),
                      str(pc_data['molecular_formula']),
                      pc_data['molecular_weight'],
                      pc_data['canonical_smiles'],
                      pc_data['isomeric_smiles'],
                      pc_data['inchi'],
                      pc_data['inchikey'],
                      pc_data['iupac_name'],
                      pc_data['xlogp'],
                      pc_data['exact_mass'],
                      pc_data['monoisotopic_mass'])

        # create the actual record in the database
        cur.execute('''INSERT INTO PC_data(cid,
                                            elements,
                                            atoms,
                                            bonds,
                                            molecular_formula,
                                            molecular_weight,
                                            canonical_smiles,
                                            isometric_smiles,
                                            inchi,
                                            inchikey,
                                            iupac_name,
                                            xlogp,
                                            exact_mass,
                                            monoisotopic_mass)
                                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', write_data)
        conn.commit()
    return True

def add_cas_data(cas_data, conn, cur):
    '''
    Adds cas data to the table if the CAS nr is not yet in the table
    '''

    # Check for existing records for this identifier in the database, a new record is only added if it does not exist yet.
    cur.execute('SELECT * from CAS_data WHERE cas_rn = ?', (cas_data['rn'],))
    result = cur.fetchall()
    
    if not result:
        # Create a tuple of the data to store in the correct order
        data = (str(cas_data['rn']),
                str(cas_data['uri']),
                str(cas_data['name']),
                str(cas_data['smile']),
                str(cas_data['canonicalSmile']),
                str(cas_data['inchi']),
                str(cas_data['inchiKey']),
                str(cas_data['molecularFormula']),
                str(cas_data['molecularMass']),
                str(cas_data['experimentalProperties']),
                str(cas_data['propertyCitations']),
                str(cas_data['synonyms']),
                str(cas_data['replacedRns']))

        # create the actual record in the database
        cur.execute('''INSERT INTO CAS_data(cas_rn,
                                            uri,
                                            name,
                                            smile,
                                            canonical_smile,
                                            inchi,
                                            inchikey,
                                            molecular_formula,
                                            molecular_weight,
                                            documented_properties,
                                            sources,
                                            synonyms,
                                            replaced_cas)
                                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''', data)
        conn.commit()
    return True

def store_data(entry_name, cas_find, pc_find, conn, cur):
    '''
    Stores all available API data in the database and makes sure the compound is referencing the correct API data records
    '''

    # Sets empty flags to be replaced by the PubChem or CAS identifier if available, if no identifier is available for either service, this remains None (thus empty).
    cas = None
    pc = None

    # Stores CAS API data if available
    if cas_find:
        add_cas_data(cas_find, conn, cur)
        cas = cas_find['rn']

    # Stores PubChem data if available
    if pc_find:
        add_pc_data(pc_find, conn, cur)
        pc = pc_find.cid

    # Add identifiers to the compound referencing the correct data records.
    # By design this is done for every compound with the same name to for efficiency
    cur.execute('UPDATE Compound_entries SET PC_data_id = ?, CAS_data_id = ? WHERE lower(compound_name) = ?', (pc, cas, entry_name.lower()))
    conn.commit()
   
def add_skiplist(entry_id):
    '''
    Adds an ID to the skiplist text file
    '''
    # create skiplist file if it doesnt exist and add the id at the end of the file
    with open("skiplist.txt",'a+') as skiplist:
        skiplist.write(entry_id+'\n')

def add_all_skiplist(entry_name, cur):
    '''
    Adds all IDs of compounds with the specified name to the skiplist
    '''
    # retrieve the id's of records with the same compound name
    cur.execute('SELECT id FROM Compound_entries WHERE lower(compound_name) = ?',(entry_name.lower(),))

    # Add those id's to the skiplist
    for i in cur.fetchall():
        add_skiplist(str(i[0]))

def check_skiplist(entry_id):
    '''
    Check if an ID is present on the skiplist
    '''
    skiplist = get_skiplist()

    # Returns True if the ID is present on the skiplist
    if entry_id in skiplist:
        return True
    
    return False

def get_skiplist():
    '''
    Returns skiplist contents as a list object
    '''

    # Opens and reads file if it exists and returns a list of ID's
    if os.path.exists("skiplist.txt"):
        with open("skiplist.txt",'r') as skiplist:
            skip = skiplist.read().splitlines()
            
    else:
        # Return an empty list if there's no skiplist file
        # the file will be created the first time an ID gets added to it
        skip = []
        
    return skip
          
def next_entry(cur, skiplist = None):
    '''
    Returns the next database record that has not yet been identified and requires identification
    '''

    # retrieves the skiplist if none was given when calling this function
    if not skiplist:
        skiplist = get_skiplist()

    # Because the skiplist has a variable length, the amount of variables to filter for should have a custom amount of '?' to use as placeholders
    skipstring = ','.join(len(skiplist)*'?')

    # SQL to retrieve one database record that has no associated CAS or PubChem ID and does not occur on the skiplist
    cur.execute('SELECT id, compound_name FROM Compound_entries WHERE (CAS_data_id IS NULL AND PC_data_id IS NULL) AND id NOT IN (%s)' % skipstring, skiplist)
    return cur.fetchone()

def run_compound(entry_id:int, entry_name:str, conn, cur):
    '''
    Automatic and manual verification of compounds using several steps:
    1. Check whether one of the compound names matches completely -> accept that match and sync the apis using InChI
    2. Prompt user for selecting the correct compound from the list
    3. Give option for manual entry of CAS or CID
    4. Give option for not being able to correctly verify a compound using any of the available means
    
    Upon verification of a compound, API data is added to the database
    '''

    # Indicate new entry and show stats
    print("===========================================")
    print(f"Current compound: {entry_name} | {compounds_to_go(cur)} left")

    # Retrieve data from both API's
    pc_data = pubchempy.get_compounds(entry_name, 'name')
    cas_data = json.loads(cas_api.search(entry_name))['results']

    # Set flags for confirmed results to False
    pc_find = False
    cas_find = False

    # CAS register Automatching
    # Loop through cas search results
    for i in cas_data:
        
        # On exact name match, use this result as verified
        if entry_name.lower() == i['name'].lower():
            print("Found CAS by name")

            # Set the found result to the result of the /details of the CAS api for this registration number
            cas_find = cas_api.details(i['rn'])
            
            # Sometimes CAS api does not have an InChI, ifso we remove undo the verification and let the user handle the situation
            if not cas_find['inchi']:
                cas_find = False
                break

            # Search the PubChem database using the InChI to match results
            pc_result = pubchempy.get_compounds(cas_find['inchi'], 'inchi')

            # InChI should only describe 1 compound, if it matches more then one this could indicate ambiguity
            if len(pc_result) == 1:
                pc_find = pc_result[0]
                print(f"PC {pc_find.to_dict()['cid']} found from inchi")

    # Store results in database on automatch of both PubChem and CAS 
    if cas_find and pc_find:
        store_data(entry_name, cas_find, pc_find, conn, cur)
        return True

    # PubChem Automatching
    # Its recommended to disable Pubchem matching due to problems described on GitHub.
    # Possibly has to do with the fact that CAS api needs the "InChI=" part removed from the InChI string.
    
    # -- COMMENT OUT TO DISABLE PubChem Automatch [START]

    # Loop through the PubChem search results for this name
    for i in pc_data:
        
        # Accept result on full IUPAC name match
        if entry_name.lower() == i.to_dict()['iupac_name'].lower():
            pc_find = i

            # Search the CAS register for the InChI.
            # this statement might be the origin of the bug; a possible untested solution might be:
            # cas_result = json.loads(cas_api.search(pc_find.to_dict()['inchi'].replace("InChI=","")))['results']
            # fixing this bug is outside the scope of the project, but could be useful for future use of the script
            cas_result = json.loads(cas_api.search(pc_find.to_dict()['inchi']))['results']

            # If a result was found throught CAS search, the chemical information is requested from the CAS API /details
            if cas_result:
                cas_find = cas_api.details(cas_result[0]['rn'])
                print(f"CAS {cas_find['rn']} found from inchi")
    
    # Store results in database on automatch of both PubChem and CAS        
    if cas_find and pc_find:
        store_data(entry_name, cas_find, pc_find, conn, cur)
        return True

    # -- COMMENT OUT TO DISABLE PubChem Automatch [END]
    
    
    # Pseudocode for the manual part of the verification process:
    # statements for:
    #   select either a cas or cid from the list of API results
    #   OR provide a manual CAS or CID
    #   OR to put on the skiplist
    # if manual entry was used or CAS/CID was selected from the list -> add the other through inchi

    # Create one list of results to print, and store information about the origin of the data
    # because both require their own method of storing due to the difference in data provided through both APIs
    total_results = []
    for result in pc_data:
        total_results.append({'result':result, 'origin':'pc'})

    for result in cas_data:
        total_results.append({'result':result, 'origin':'cas'})

    # Information for the user about how to progress
    print("Type index number and press enter to confirm")

    # Show the user the third-party identifier and the name provided through the API
    for index, result in enumerate(total_results):
        if result['origin'] == 'pc':
            ident = result['result'].cid
            name = result['result'].iupac_name
        else:
            ident = result['result']['rn']
            name = result['result']['name']
        print(f'{index}:\t{ident}\t-\t{name}')

    while True:
        # Request the user to choose an index of a result or one of the three other options (skiplist, manual pubchem entry or manual CAS entry) 
        print("Or press 's' to add to skiplist, 'mp' for manual PubChem CID, 'mc' for manual CAS number")
        choice = input('>')

        # == SELECTED LISTED RESULT == 
        # If the choice is a numerical, it indicates the user chose one of the results.
        if choice.isdigit():
            # user chose one of the proposed entries
            chosen = total_results[int(choice)]
            
            # PubChem chosen
            if chosen['origin'] == 'pc':
                # Store the chosen PubChem data
                pc_find = chosen['result']

                # Search the CAS register using the InChIKey. It seems that I ran into problems using InChI
                # presumably due to the the InChI= part still being in front of the InChI which is not included in the InChIKey.
                # This did work well however
                cas_result = json.loads(cas_api.search(pc_find.to_dict()['inchikey']))['results']
                if cas_result:
                    cas_find = cas_api.details(cas_result[0]['rn'])
                    print(f"CAS {cas_find['rn']} found from inchi")

            # CAS chosen
            else:
                # Store the chosen CAS data
                cas_find = cas_api.details(chosen['result']['rn'])

                # search PubChem by InChI
                pc_result = pubchempy.get_compounds(cas_find['inchi'], 'inchi')
                if len(pc_result) == 1:
                    pc_find = pc_result[0]
                    print(f"PC {pc_find.to_dict()['cid']} found from inchi")

            # Store whatever APIs returned data for the chosen compound into the Database and end the function
            store_data(entry_name, cas_find, pc_find, conn, cur)
            return True

        # == SELECTED SKIPLIST ==
        elif choice.lower() == "s":
            # add to skiplist
            add_all_skiplist(entry_name, cur)
            return False

        # == SELECTED MANUAL PUBCHEM ENTRY ==
        elif choice.lower() == "mp":

            # Prompt user to enter a CID and make sure it's correctly formatted
            while True:
                print('Please paste full PubChem CID here and press enter to confirm')
                cid = input('>').strip()
                if cid.isdigit():
                    break

            # Search PubChem for the CID
            pc_result = pubchempy.get_compounds(cid, 'cid')

            # Accept the result if this indeed only resulted in one option
            if len(pc_result) == 1:
                pc_find = pc_result[0]

                # use the InChI to search the CAS register --> this due to it being an InChI including "InChI=" might again result in no CAS results
                # that's why it was opted to not use the manual PubChem unless no suitable CAS number could be found manually
                cas_result = json.loads(cas_api.search(pc_find.to_dict()['inchi']))['results']
                if cas_result:
                    cas_find = cas_api.details(cas_result[0]['rn'])
                    print(f"CAS {cas_find['rn']} found from inchi")

            # Store the API data in the database for whatever API returned data
            store_data(entry_name, cas_find, pc_find, conn, cur)
            return True
                
        # == SELECTED MANUAL CAS ENTRY ==
        elif choice.lower() == "mc":

            # Prompt user for the CAS number
            while True:
                print('Please paste full CAS registration number here and press enter to confirm')
                cas = input('>').strip()

                # Requirements for CAS numbers are harder than this, but the goal here was to prevent accidental input
                if len(cas) > 4:
                    break

            # Use the user input to request API for the /details page
            cas_result = cas_api.details(cas)

            # If the CAS details are available store the CAS result
            if cas_result:
                cas_find = cas_result

                # Use the InChI to search the cas register
                pc_result = pubchempy.get_compounds(cas_find['inchi'], 'inchi')

                # Only one compound should match this InChI
                if len(pc_result) == 1:
                    pc_find = pc_result[0]
                    print(f"PC {pc_find.to_dict()['cid']} found from inchi")

            # Store the API data of whatever API returned data
            store_data(entry_name, cas_find, pc_find, conn, cur)
            return True
        

def compounds_to_go(cur, skiplist=None):
    '''
    Prints amount of entries that require data and are not present on the skiplist
    '''
    # This functions very similar to retrieving the next entry
    # if no skiplist was provided, retrieve the current skiplist
    if not skiplist:
        skiplist = get_skiplist()

    # create a variable length string of '?' to allow for the amount of variables in the SQL statements as there are IDs in the skiplist
    skipstring = ','.join(len(skiplist)*'?')

    # Retrieve the count of distinct compound names that have no associated CAS or PubChem data and are not present on the skiplist
    cur.execute('SELECT count(DISTINCT lower(compound_name)) FROM Compound_entries WHERE (CAS_data_id IS NULL AND PC_data_id IS NULL) AND id NOT IN (%s)' % skipstring, skiplist)
    return cur.fetchone()[0]


if __name__ == "__main__":
    
    # Setup the database connection
    conn,cur = db()
    
    while True:
        # Retrieve the next entry
        entry_id, entry_name = next_entry(cur)

        # if no more entries are available end the loop
        if not entry_name:
            print("No more entries.")
            break

        # if an entry is available procss the entry
        run_compound(entry_id, entry_name, conn, cur)
