#!/usr/bin/python3

import urllib.parse
import urllib.request

def uniprot2ids(uniprotid,to,query):
  params = {
  'from': uniprotid,
  'to': to,
  'format': 'tab',
  'query': query
  }
  data = urllib.parse.urlencode(params)
  data = data.encode('utf-8')
  url = 'https://www.uniprot.org/uploadlists/'
  req = urllib.request.Request(url, data)
  with urllib.request.urlopen(req) as f:
     response = f.read()
  r = response.decode('utf-8')
  return(r)
