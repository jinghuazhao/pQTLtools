#!/usr/bin/python3

import requests, sys

requestURL = 'https://www.ebi.ac.uk/proteins/api/coordinates?offset=0&size=100&gene=FGFR2'
r = requests.get(requestURL, headers={ "Accept" : "application/json"})
if not r.ok:
   r.raise_for_status()
   sys.exit()
responseBody = r.text
print(responseBody)
