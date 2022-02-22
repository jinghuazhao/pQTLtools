#!/usr/bin/bash

import requests, sys
requestURL = "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession=P21802"
r = requests.get(requestURL, headers={ "Accept" : "application/xml"})
if not r.ok:
   r.raise_for_status()
   sys.exit()
responseBody = r.text
print(responseBody)
