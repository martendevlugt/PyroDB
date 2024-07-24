import re
import sqlite3
from identifier import db

def result_to_float(result):
    ''' Returns a float value if a given value is not empty, if it is empty it returns a float zero'''

    # empty values are normally given as multiple dashes (always the same) these are turned into a float zero, otherwise a float of the extracted text is returned
    if result == "------":
        return 0.0
    else:
        return float(result)

def extract_epi_summary(raw_data):
    '''
    Extract Summary from EPI suite output file
    '''

    # Isolate the summary from raw_data
    pattern = re.compile(r"-{30} EPI SUMMARY.*-{26}(.*?    Melting Point \(deg C\)  :   .*?\n)", re.DOTALL)
    data = pattern.search(raw_data).group(0)

    # Setting up patterns
    solubility_pattern = re.compile(r"    Water Solubility \(mg/L\):   (.*?)\n", re.DOTALL)
    vapor_pattern = re.compile(r"    Vapor Pressure \(mm Hg\) :   (.*?)\n", re.DOTALL)
    henry_pattern = re.compile(r"    Henry LC \(atm-m3/mole\) :   (.*?)\n", re.DOTALL)
    kow_pattern = re.compile(r"    Log Kow \(octanol-water\):   (.*?)\n", re.DOTALL)
    boiling_pattern = re.compile(r"    Boiling Point \(deg C\)  :   (.*?)\n", re.DOTALL)
    melting_pattern = re.compile(r"    Melting Point \(deg C\)  :   (.*?)\n", re.DOTALL)

    # extracting data using patterns
    processed = {}
    processed['solubility'] = result_to_float(solubility_pattern.search(raw_data).group(1))
    processed['vapor'] = result_to_float(vapor_pattern.search(raw_data).group(1))
    processed['henry'] = result_to_float(henry_pattern.search(raw_data).group(1))
    processed['kow'] = result_to_float(kow_pattern.search(raw_data).group(1))
    processed['boiling'] = result_to_float(boiling_pattern.search(raw_data).group(1))
    processed['melting'] = result_to_float(melting_pattern.search(raw_data).group(1))

    return processed

def split_compounds(raw_data):
    '''
    Splits EPI suite output data into the model outputs per compound.
    Returns a list of individual compound outputs.
    '''
    compounds = raw_data.split('''


========================


''')
    return compounds

def extract_ecosar(raw_data):
    '''
    Extracts Ecosar model data from EPI suite output file
    '''

    # Retrieve ECOSAR version info
    pattern = re.compile(r"ECOSAR Version(.*?)Note:.*?\n", re.DOTALL)
    data = pattern.search(raw_data).group(0)

    # Base parameter patterns
    smiles_pattern = re.compile(r"SMILES : (.*?)\n", re.DOTALL)
    chem_pattern = re.compile(r"CHEM   : (.*?)\n", re.DOTALL)
    cas_pattern = re.compile(r"CAS Num: (.*?)\n", re.DOTALL)
    chemid_pattern = re.compile(r"ChemID1: (.*?)\n", re.DOTALL)
    molfor_pattern = re.compile(r"MOL FOR: (.*?)\n", re.DOTALL)
    molwt_pattern = re.compile(r"MOL WT : (.*?)\n", re.DOTALL)

    # Used parameters patterns
    modelparams_pattern = re.compile(r"Values used to Generate(.*)\n\n\n", re.DOTALL)
    modelparams = modelparams_pattern.search(data).group(0)
    kow_pattern = re.compile(r"Log Kow: (.*?) .*?\((.*?)\)\n", re.DOTALL)
    solubility_pattern = re.compile(r"Wat Sol: (.*?) .*?\((.*?)\)\n", re.DOTALL)

    # Extracting data using patterns and storing the processed data in a dict
    processed = {}
    processed['SMILES'] = smiles_pattern.search(data).group(1)
    processed['CHEM'] = chem_pattern.search(data).group(1)
    processed['CAS'] = cas_pattern.search(data).group(1)
    processed['ChemID'] = chemid_pattern.search(data).group(1)
    processed['mol_for'] = molfor_pattern.search(data).group(1)
    processed['mol_weight'] = result_to_float(molwt_pattern.search(data).group(1))
    
    processed['kow'] = result_to_float(kow_pattern.search(modelparams).group(1))
    processed['kow_model'] = kow_pattern.search(modelparams).group(2)
    processed['solubility'] = result_to_float(solubility_pattern.search(modelparams).group(1))
    processed['solubility_model'] = solubility_pattern.search(modelparams).group(2)

    # Extracting data for the simulated animal tests
    tests_processed = {}

    # tests are isolated from the output, the part starts with multiple '=' and always ends in a 'Note'.
    # can have a variable amount of lines in between depending on the amount of tests that were simulated.
    tests_pattern = re.compile(r"=====.*?\n(.*) Note", re.DOTALL)
    tests = tests_pattern.search(data).group(0)

    # Each test result is on a new line
    tests = tests.split('\n')
    
    tests_list = []
    for test in tests:
        try:
            # The tests are listed in a table format, cells in this table are devided using multiple spaces
            splittest = test.split('   ')

            # Cleaning up of extra spaces and other information
            parts_list = []
            for part in splittest:
                cleaned_part = part.replace(":","")
                cleaned_part = cleaned_part.strip()
                if len(cleaned_part) > 0:
                    parts_list.append(cleaned_part)

            # Required information is extracted from the cleaned table
            test_results = {}
            test_results['ecosar_class'] = parts_list[0]
            test_results['organism'] = parts_list[1]
            test_results['duration'] = parts_list[2]
            test_results['endpoint'] = parts_list[3]

            # Sometimes asterixes are used to provide a note about solubility being to low to express the predicted effects in real tests.
            # This information is extracted, but not used further in the PyroDB project.
            if '*' in parts_list[4]:
                test_results['solubility_error'] = True

            # The asterix is cleaned off when it was present after flagging its existence
            test_results['predicted_conc'] = result_to_float(parts_list[4].replace(" *",""))

            # results are stored
            tests_list.append(test_results)
        except:
            pass

    # store test results list in the output dict
    processed['tests'] = tests_list
    return processed

def extract_biowin(raw_data):
    '''
    Extract BioWin data from EPI suite output data returns results for each of the 7 models and the ready biodegradation prediction
    '''

    # Isolate BioWin data from raw_data
    data_pattern = re.compile(r"BIOWIN .*? Program Results:(.*?)\n\n\n\n",re.DOTALL)
    data = data_pattern.search(raw_data).group(0)

    # Set up patterns for extracting the conclusion of each of the models included in BioWin
    biowin1_pattern = re.compile(r"Biowin1 \(Linear Model Prediction\)    :  (.*?)\n",re.DOTALL)
    biowin2_pattern = re.compile(r"Biowin2 \(Non-Linear Model Prediction\):  (.*?)\n",re.DOTALL)
    biowin3_pattern = re.compile(r"Biowin3 \(Ultimate Biodegradation Timeframe\):  (.*?)\n",re.DOTALL)
    biowin4_pattern = re.compile(r"Biowin4 \(Primary  Biodegradation Timeframe\):  (.*?)\n",re.DOTALL)
    biowin5_pattern = re.compile(r"Biowin5 \(MITI Linear Model Prediction\)    :  (.*?)\n",re.DOTALL)
    biowin6_pattern = re.compile(r"Biowin6 \(MITI Non-Linear Model Prediction\):  (.*?)\n",re.DOTALL)
    biowin7_pattern = re.compile(r"Biowin7 \(Anaerobic Model Prediction\):  (.*?)\n",re.DOTALL)
    ready_pattern = re.compile(r"Ready Biodegradability Prediction:  (.*?)\n",re.DOTALL)

    # Use patterns to retrieve model conslusions, store these with model numbers
    results = {}
    results['1'] = biowin1_pattern.search(data).group(1)
    results['2'] = biowin2_pattern.search(data).group(1)
    results['3'] = biowin3_pattern.search(data).group(1)
    results['4'] = biowin4_pattern.search(data).group(1)
    results['5'] = biowin5_pattern.search(data).group(1)
    results['6'] = biowin6_pattern.search(data).group(1)
    results['7'] = biowin7_pattern.search(data).group(1)
    results['ready'] = ready_pattern.search(data).group(1)

    # The model conclusions are based on numerical model outputs, these are stored for each model as [model_number][_value]
    # Outputs look as follows, the desired value is (if split at the '|' character) at index 3:
    # RESULT   |   Biowin7 (Anaerobic Linear Biodeg Prob)   |         |  0.7555
    # Each of the models contains such output in order of appearance with the exclusion of ready biodegradability.
    # I'm sorry for the horrible name of the pattern variable name
    testtest = re.compile(r"(RESULT.*?)\n",re.DOTALL)
    for number, i in enumerate(testtest.findall(data)):
        results[f"{number+1}_value"] = result_to_float(i.split('|')[3].strip())

    return results

def extract_bcfbaf(raw_input):
    '''
    Extracts BCFBAF model values from EPI suite output
    '''

    # BCF and BAF patterns
    # All useable values are captured in 2 lines, both the log and and the non log version
    bcf_pattern = re.compile("Log BCF \(regression-based estimate\):(.*?)\(BCF = (.*?) ", re.DOTALL)
    baf_pattern = re.compile("Log BAF \(Arnot-Gobas upper trophic\):(.*?)\(BAF = (.*?) ", re.DOTALL)

    # Extracting data using defined patterns
    result = {}
    result['log_bcf'] = result_to_float(bcf_pattern.search(raw_input).group(1).strip())
    result['bcf'] = result_to_float(bcf_pattern.search(raw_input).group(2).strip())
    result['log_baf'] = result_to_float(baf_pattern.search(raw_input).group(1).strip())
    result['baf'] = result_to_float(baf_pattern.search(raw_input).group(2).strip())
    
    return result


def assessment(test_results):
    '''
    Runs a full assessment on data printing the result and returning only persistence, bioaccumulativity and toxicity
    '''
    persistence = get_persistence(test_results['biowin'])
    bioaccumulativity = get_bioaccumulativity(test_results['bcfbaf'])
    toxicity = get_toxicity(test_results['ecosar'])
    ready_soluble = get_solubility(test_results['ecosar'])

    # I noticed a difference when running several variations of smiles for the same compound, and discovered that sometimes EPI suite uses CAS supplied data to get better results (like kow etc)
    # I made a flag to keep track of this fact
    dbstat = "stor" if test_results["base_info"]["using_db"] else "calc"

    print(f'Screening result ({test_results["base_info"]["id"]}) {dbstat}: {persistence}\t{bioaccumulativity}\t{toxicity}\t{ready_soluble}\t{test_results["ecosar"]["SMILES"]}')
    return (persistence, bioaccumulativity, toxicity)

def get_solubility(ecosar):
    '''
    Assesses solubility using the Kow used to run the ECOSAR model
    A compound is deemed soluble be low the log(Kow) of octanol (3)
    '''
    ready_soluble = 'S' if ecosar['kow'] < 3 else ''
    return ready_soluble

def get_persistence(biowin):
    '''
    Assessment for persistence using biowin model data
    Chosen parameters:
    page 60: https://echa.europa.eu/documents/10162/13632/information_requirements_r11_en.pdf/a8cce23f-a65a-46d2-ac68-92fee1f9e54f
    biowin 2 probability < 0.5 means biodegrades fast
    biowin 3 value < 2.25 (to 2.75) chosen to go for high value, because these would normally need additional screening
    biowin 6 probability < 0.5 means biodegrades fast
    '''
    persistence = ""
    if biowin['2_value'] < 0.5 and biowin['3_value'] < 2.75:
        persitence = "P/vP"
    if biowin['6_value'] < 0.5 and biowin['3_value'] < 2.75:
        persistence = "P/vP"
    return persistence

def get_bioaccumulativity(bcfbaf):
    '''
    Assessment for bioaccumulativity based on values from the BCFBAF model
    B/VB on basis of BCFBAF model, page 20 of echa manual
    BCF of 2000 implies bioaccumulative, over 4000 implies very bioaccumulative
    '''
    bioaccumulativity = ''
    if bcfbaf['bcf'] > 2000:
        bioaccumulativity = 'B'
    if bcfbaf['bcf'] > 4000:
        bioaccumulativity = 'vB'
    return bioaccumulativity

def get_toxicity(ecosar):
    '''
    Assesses toxicity using ECOSAR model outputs
    page 134 defines that toxicity occurs when either LC50 or EC50 is less than 0.01 mg/l
    '''
    max_concentration = 10000 # set unreasonably high otherwise the condition will never trigger

    for test in ecosar['tests']:
        if (test['organism'] in ['Fish', 'Daphnid', 'Green Algae'] and
            test['endpoint'] in ['LC50','EC50'] and
            test['predicted_conc'] <= max_concentration):
            max_concentration = test['predicted_conc']
    
    toxicity = 'T' if max_concentration <= 0.01 else ''
    return toxicity

def extract_base(raw_data):
    '''
    Extracts some base information like the user given ID of a compound and whether CAS data was used during the run
    '''

    # Setting up patterns
    entryid_pattern = re.compile("CHEM   : (.*?)\n", re.DOTALL)
    usingdb_pattern = re.compile("CAS Num  :  ")

    # Extracting data from patterns
    entry_id = entryid_pattern.search(raw_data).group(1)
    usingdb_result = usingdb_pattern.search(raw_data)

    # Turning the result from usingdb into a True/False flag, a cas number is only presented through EPI suite if it uses data from the internal database
    usingdb = True if usingdb_result else False
    
    return {'id':entry_id, 'using_db':usingdb}

def get_inchi(ident:str, ident_file = "translation.txt"):
    '''
    Translates the user specified ID back into the InChI of the compound
    '''
    # Open and read the translation file
    with open(ident_file) as file:
        translation = file.read().splitlines()

    # Find the correct ID, and return associated InChI
    for i in translation:
        i = i.split('\t')
        if i[0] == ident:
            return i[1]

def store_result(inchi, result, conn, cur):
    '''
    Stores assessment results and data used to get the result into the dataset.db file
    '''

    # Turns the True/False into 1/0 respectively
    using_stored = 1 if result['base_info']['using_db'] else 0

    # Recalculates the lowest concentration which was used to base toxicity on (suboptimal, but gets the job done)
    max_concentration = 10000 # set unreasonably high otherwise the condition will never trigger
    for test in result['ecosar']['tests']:
        if (test['organism'] in ['Fish', 'Daphnid', 'Green Algae'] and
            test['endpoint'] in ['LC50','EC50'] and
            test['predicted_conc'] <= max_concentration):
            max_concentration = test['predicted_conc']

    # Creates a tuple of the data in the correct order
    data = (inchi,
            get_persistence(result['biowin']),
            get_bioaccumulativity(result['bcfbaf']),
            get_toxicity(result['ecosar']),
            get_solubility(result['ecosar']),
            using_stored,
            result['bcfbaf']['bcf'],
            max_concentration,
            result['biowin']['2_value'],
            result['biowin']['3_value'],
            result['biowin']['6_value'],
            result['ecosar']['kow'])

    # Stores the assessment results and important values in the database
    cur.execute("""INSERT INTO Ecotoxicity(inchi,P,B,T,S,using_stored,BCFBAF,ECOSAR,BIOWIN2,BIOWIN3,BIOWIN6,solubility) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""", data)
    conn.commit()

def main(infile):
    '''
    main process running for processing EPI suite data
    '''

    # setup database connection
    conn, cur = db()

    # get the full output file as one large string
    with open(infile,'r') as f:
        raw_data = f.read()

    # split the full output string into a list with an entry for each individual run
    compound_tests = split_compounds(raw_data)

    last_id = '0'
    found_stored = False
    last_result = []

    for test in compound_tests:
        # Set up a dictionary to store model results
        test_results = {}

        # Retrieve data for each model
        test_results['base_info'] = extract_base(test)
        test_results['epi_summary'] = extract_epi_summary(test)
        test_results['ecosar'] = extract_ecosar(test)
        test_results['biowin'] = extract_biowin(test)
        test_results['bcfbaf'] = extract_bcfbaf(test)

        # Retrieve the ID corresponding to the EPI SMILES batch file
        current_id = test_results['base_info']['id'].split("_")[0]

        # pseudocode explaining how the data to store into the database was chosen:
        # for each new id:
        #   if this id data used stored info, store this info an stop processing this id by setting found flag true
        #   if we found no stored and we start processing a new id now store the old one
        #   if we're processing a new id and the found store flag is still true, reset it to false

        # If no test using data from the internal EPI suite database was found, we store the last result
        if not found_stored and current_id != last_id:
            
            # Show assessment data, uncomment line below to enable output on screen
            assessment(last_result)

            found_stored = False
            
            # Store result to the database, uncomment line below to stop storing in the database
##            store_result(get_inchi(last_id), last_result, conn, cur)

        # Reset the flag that indicates a test using data from the EPI internal database was found when processing a new test
        elif found_stored and current_id != last_id:
            
            found_stored = False

        # When a test using data from the internal EPI suite database is found it is stored and the rest of the tests for the same compound id are ignored
        if test_results['base_info']['using_db'] and not found_stored:
            # Uncomment the line below to show assessment data on screen
            assessment(test_results)

            # Uncomment below to prevent storing data in the compound entry database
##            store_result(get_inchi(current_id), test_results, conn, cur)

            # Set flag to ignore other tests on the same compound
            found_stored = True

        # store the data of this test for one cycle for when we need to store the previous test
        last_result = test_results
        last_id = current_id
##        assessment_result = assessment(test_results)

if __name__ == "__main__":
    # specifies the EPI suite full output file to use
    infile = "new_results.OUT"

    # run the process to automatically add the EPI suite results to the compound entry database
    main(infile)
