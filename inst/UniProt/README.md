# Proteins API

Source [https://www.ebi.ac.uk/proteins/api](https://www.ebi.ac.uk/proteins/api)

```bash
# accession
curl -X GET --header 'Accept:application/json' https://www.ebi.ac.uk/proteins/api/proteomics?offset=0&size=-1&accession=O14788
curl -X GET --header 'Accept:application/json' https://www.ebi.ac.uk/proteins/api/proteomics/O14788

# gene
curl -X GET --header 'Accept:application/json' https://www.ebi.ac.uk/proteins/api/coordinates?gene=FGFR2

# interaction
curl -X GET --header 'Accept:application/json' https://www.ebi.ac.uk/proteins/api/proteins/interaction/O14625
```
