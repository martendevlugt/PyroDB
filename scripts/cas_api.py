import requests
import argparse
import json

def search(query):
    '''Search cas library by Cas Nr. Smiles, InChl(without prefix), InChlKey or name'''
    base_url = "https://commonchemistry.cas.org/api/search"
    full_url = f'{base_url}?q={query}'
    result = requests.get(full_url)
    if result.status_code == 200:
        return result.text
    print(result.status_code)

def details(query):
    ''' Returns details of a chemcal by cas number '''
    base_url = "https://commonchemistry.cas.org/api/detail"
    full_url = f'{base_url}?cas_rn={query}'
    result = requests.get(full_url)
    result.encoding = 'UTF-8'
    if result.status_code == 200:
        return json.loads(result.text)
    print(result.status_code)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--search', type=str, help="Name or cas nr. of a chemical to search for")
    args = parser.parse_args()

    query = args.search

    if not args.search:
        print('Please enter the name of the chemical you would like to search')
        query = input('>')

    search_result = json.loads(search(query))['results']
   
    for index, result in enumerate(search_result):
        print(f"{index+1:4d}: {result['name']:20s} - {result['rn']}")

    print("Type the index of the result you want to retrieve the details of or press q or enter to exit")
    detail_index = input(">")
    if detail_index == None or not detail_index.isdigit()or int(detail_index) > len(search_result):
        print("Invalid input or chose to exit, goodbye!")
        exit(1)

    detailed_result = details(search_result[int(detail_index)-1]['rn'])
    for i in detailed_result:
        if i != "image":
            print(f"{i}: {detailed_result[i]}")
