# PyroDB
This repository holds the code written as part of my thesis on the hazard profiles of pyrolysis plant watewater.
For this project, a database to collect compounds names from pyrolysis products was designed, together with tools to match the compound names to the CAS register and PubChem database.
Data from the CAS register and PubChem database was used to screen the compounds for possibly hazardous compounds using EPI suite.

## [Dataset](dataset.db)
An empty "dataset.db" is included which includes the correct tables and collumns. In the thesis compound names originating from research papers containing py-GC/MS analytical data were used as a starting point.
Chemical data was then retrieved using the "auto_identifier.py" and "identifier.py" scripts. Then the chemical data was used with the "epi_input.py" script to prepare for the usage of EPI suite to generate data for the screening level ecotoxicological assessment which was processed with the "epi_processor.py" script.

The database was created using [DB-browser for SQLite version 3.12.2](https://sqlitebrowser.org/dl/)

## Included scripts
### [cas_api.py](scripts/cas_api.py)
Wrapper around the CAS register API, used in the compound identification script but can be ran as a standalone script to retrieve information on a single compound.

The "cas_api.py" script has two functions:
* search
  * Takes a search query (name, InChI key, etc.) and returns a list of results
* details
  * Takes a CAS registration number (can be obtained from search) and returns information on the compound
 
### [compound_loader.py](scripts/compound_loader.py)
Allows insertion of compound names and experiment ID's from a text file into the "dataset.db" database.
Reason for this script is to convert a compound name list that was created in Excel quickly to the database via a .txt file. This was needed due to the first version being created in Excel.
Text file input should be structured as:
`[compound_name][tab][experiment_id]`

### [identifier.py](scripts/identifier.py) & [auto_identifier.py](scripts/auto_identifier.py)
Both these scripts share a part of their functionality. Which is why they share this chapter.

The identifier script uses the PubChem and CAS api's to identify and retrieve information for each distinc compound name in the "dataset.db" database. To achieve this, the [PubChemPy library](https://pubchempy.readthedocs.io) and [CAS api](scripts/cas_api.py) are used.
When finding an exact name match on either one of the API's the chemical data is added to the "dataset.db" database tables for the respective service and the compounds entry recieves the PubChem CID or CAS registration number as a reference to the retrieved data.
If no exact match was found, when using the "auto_identifier.py" script the compound is skipped. When using the "identifier.py" however, the found results are listed allowing the user to choose the correct compound. The user should now do their own research and verify that either one of the listed compounds is indeed correct; or supply a different PubChem CID or CAS registration number. Due to an inconsistency with the CAS search API please read the warning at the end of this chapter very carefully!
When a match is made either automatically or manually added, the InChI is retrieved from the chosen result and used to query the API of the service that was not chosen. This should whenever available retrieve information for the exact same compound from the other API without further user intervention.

The reason for the existense of the automatic script is to allow a pass over the full list of compounds automatically before running the manual script. This prevents cases where the user has to sit and wait (due to the time required for loading each result) idle until the script requires manual input. Making the time spent identifying compounds manually as efficient as possible with as little waiting time as possible.
Depending on the amount of compounds, the identification process can take up to several hours (1244 compounds can take up between 6 and 8 hours).

Whenever no match is found by the automatically matching or manual entry, the user can skip this compound by adding it's id to a list of database id's to skip (listed in the "skiplist.txt" file) by entering 's' when selecting a compound. 

> [!WARNING]
> The CAS api does not seem to support searching with an InChI as a query.
> 
> This leads to the script not being able to retrieve data form the CAS register when using a result provided by PubChem (either from manual entry of CID ("mp") or automatically matching on an exact name match).
> 
> To prevent this, exact name matching using PubChem should be disabled and whenever available CAS registration numbers (recognizable by the two '-' in the number) should be chosen from the results or manually added using the "mc" functionality.

### [episuite_input.py](scripts/episuite_input.py)
Used to create a file that can be used as input for [EPI suite](https://www.epa.gov/tsca-screening-tools/download-epi-suitetm-estimation-program-interface-v411)
Automatically creates 2 files:
* epi_input.txt
  * File that can be used as a Batch Smiles (F5) input file in EPI suite. This file will contain an ID and every type of SMILE available for a compound
  * Format: `[SMILE][space][ID][_xy]` ('x' denotes the the origin of the smile: 'p' for PubChem and 'c' for CAS. 'y' denotes the type of smile: 's' for smile, 'i' for isomeric and 'c' for canonical)
* translation.txt
  * Translates the ID to an InChI corresponding to a compound from the "dataset.db"

The combination of these files allows tracking which results belong to which exact compound by allowing the result to be matched to a compound InChI.
The reasoning behind this is that when running a SMILE through EPI suite, EPI suite tends to reformat the smile. This makes it harder to match the result of EPI suite to the correct compound, whilst due to having the InChI's we have a better way to identify a compound.
The actual InChI string however is quite long and contains characters that could cause unknown errors in EPI suite. To prevent running into bugs it was decided that abstracting the InChI's to a smaller amount of numbers and letters would probably be wise.
Another choice that was made was to run all the variations of smiles available through EPI suite. EPI suite seems to use an internal database that has many SMILES stored with actual real world data for these compounds and prefers to use this real world data whenever it recognizes a SMILE. The matching however seems to be done on the literal untransformed SMILE used as input.
To make sure that if any real world data is present about a compound we actually use this, we run all available smiles hoping at least one of them hits. However this increases the amount of compounds that are ran through EPI suite by at least 2 to 4 times (both CAS and PubChem supply up to 2 variations of SMILES).

### [epi_processor.py](scripts/epi_processor.py)
EPI suite should be ran with output set to FULL in batch mode.

EPI suite creates a ".OUT" file containing the model output for each of the models included within EPI suite in one large file.
The "epi_processor.py" script takes this output files and used the RE (Regular Expressions) library to extract the relevant data from the model outputs. Then it runs the relevant data through a set of rules for screening based on ECHA "Guidance on Information Requirements
and Chemical Safety Assessment" [Chapter R.11: PBT/vPvB assessment](https://www.echa.europa.eu/documents/10162/17224/information_requirements_r11_en.pdf)

This script will either print assessment results or store the results directly into the "ecotoxicology" table in the "dataset.db" database or do both. Both these functions can be enabled/disabled by commenting out their respective lines in the "main()" function (the comments in the code will tell you which lines this refers to)
Note: This script returns one assessment result per ID (see the chapter "episuite_input.py"). If it is confirmed that one of the results uses real world data, this is the result that will be used; otherwise the results should all be the same and it will use the last result in order of occurence in the ".OUT" file.
