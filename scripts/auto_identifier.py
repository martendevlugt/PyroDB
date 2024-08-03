import identifier #uses the identifier.py script for db function; could have used more common functions
import pubchempy
import cas_api
import json

def next_entry(cur, skip):
    '''
    Returns the next database entry that has not yet been identified and requires identification
    '''

    # Because the skiplist has a variable length, the amount of variables to filter for should have a custom amount of '?' to use as placeholders
    skipstring = ','.join(len(skip)*'?')

    # SQL to retrieve one database record that has no associated CAS or PubChem ID and does not occur on the skiplist
    cur.execute('SELECT id, compound_name FROM Compound_entries WHERE (CAS_data_id IS NULL AND PC_data_id IS NULL) AND id NOT IN (%s)' % skipstring, skip)
    return cur.fetchone()

def add_pc_data(pc_data, conn, cur):
    '''
    Add pubchem data if this cid is not yet in the dataset
    '''

    # Creates dict version from the PubChemPy Compound object
    pc_data = pc_data.to_dict()

    # Check for existing records for this identifier in the database, a new record is only added if it does not exist yet.
    cur.execute('SELECT * FROM PC_data WHERE cid = ?', (pc_data['cid'],))
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
                                        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)''',data)
        conn.commit()
    return True

def store_data(entry_name, cas_find, pc_find, conn, cur):
    '''
    Stores the retrieve data from any of the API's in the database
    '''
    cas = None
    pc = None
    if cas_find:
        add_cas_data(cas_find, conn, cur)
        cas = cas_find['rn']
    if pc_find:
        add_pc_data(pc_find, conn, cur)
        pc = pc_find.cid

    cur.execute('UPDATE Compound_entries SET PC_data_id = ?, CAS_data_id = ? WHERE lower(compound_name) = ?', (pc, cas, entry_name.lower()))
    conn.commit()
        
def run_compound(entry_name, conn, cur):
    '''
    Only automatic verification of compounds
    Checks PubChem and CAS API to match name or IUPAC name (whichever available), if an exact match is found this is considered the correct compound
    '''
    # Retrieve first data due to incosistency in earlier version, should only affect retrieval of data.
    entry = entry_name

    # Retrieve data from both API's
    pc_data = pubchempy.get_compounds(entry, 'name')
    cas_data = json.loads(cas_api.search(entry))['results']

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

            # Search the PubChem database using the InChI to match results
            pc_result = pubchempy.get_compounds(cas_find['inchi'], 'inchi')
            
            # InChI should only describe 1 compound, if it matches more then one this could indicate ambiguity
            if len(pc_result) == 1:
                pc_find = pc_result[0]
                print(f"PC {pc_find.to_dict()['cid']} found from inchi")

    # Store results in database on automatch of both PubChem and CAS 
    if cas_find and pc_find:
        store_data(entry_name, cas_find, pc_find, conn, cur)
        return "full"

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

    
                
    if cas_find and pc_find:
        store_data(entry_name, cas_find, pc_find, conn, cur)
        return "full"
    
    # -- COMMENT OUT TO DISABLE PubChem Automatch [END]

    # Report back a partial find as partial
    if cas_find or pc_find:
        store_data(entry_name, cas_find, pc_find, conn, cur)
        return "partial"

    # If no matches where found we report back a failure
    return "failure"

if __name__ == "__main__":
    # Setup database connection
    conn, cur = identifier.db()

    # Create local skiplist
    skiplist = []

    # Keep stats to show progress
    partial_success = 0 # One API verified
    full_success = 0 # Both APIs verified
    failures = 0 # Neither API verified
    
    to_go = identifier.compounds_to_go(cur)

    # Run until all compounds are processed
    while True:
        # Get next compound
        entry_id, entry_name = next_entry(cur, skiplist)

        # End if there are no more compounds
        if not entry_name:
            break

        # Show stats
        print(f"{entry_name} tg={identifier.compounds_to_go(cur)}, {partial_success=}, {full_success=}, {failures=}")

        # Run identification steps
        status = run_compound(entry_name, conn, cur)

        # Update stats
        match status:
            case "full":
                full_success += 1
            case "partial":
                partial_success += 1
            case _:
                failures += 1
                skiplist.append(entry_id)
                
        # Update compounds to go
        to_go = identifier.compounds_to_go(cur)
