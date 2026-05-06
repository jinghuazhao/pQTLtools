# Protein-altering variants

Variants which directly affect the coding sequence or the protein
produced.

## Usage

``` r
data(protein_altering_variants)
```

## Format

A vector of consequences.

## Examples

``` r
protein_altering_variants
#>  [1] "coding_sequence_variant"                                      
#>  [2] "frameshift_variant"                                           
#>  [3] "frameshift_variant,NMD_transcript_variant"                    
#>  [4] "frameshift_variant,splice_donor_region_variant"               
#>  [5] "frameshift_variant,splice_region_variant"                     
#>  [6] "frameshift_variant,splice_region_variant,intron_variant"      
#>  [7] "frameshift_variant,start_lost,start_retained_variant"         
#>  [8] "frameshift_variant,stop_lost"                                 
#>  [9] "frameshift_variant,stop_retained_variant"                     
#> [10] "incomplete_terminal_codon_variant,coding_sequence_variant"    
#> [11] "inframe_insertion"                                            
#> [12] "inframe_insertion,NMD_transcript_variant"                     
#> [13] "inframe_insertion,splice_region_variant"                      
#> [14] "inframe_insertion,stop_retained_variant"                      
#> [15] "missense_variant"                                             
#> [16] "missense_variant,NMD_transcript_variant"                      
#> [17] "missense_variant,splice_donor_region_variant"                 
#> [18] "missense_variant,splice_region_variant"                       
#> [19] "missense_variant,splice_region_variant,NMD_transcript_variant"
#> [20] "missense_variant,stop_retained_variant"                       
#> [21] "protein_altering_variant"                                     
#> [22] "stop_gained"                                                  
#> [23] "stop_gained,frameshift_variant"                               
#> [24] "stop_gained,frameshift_variant,splice_region_variant"         
#> [25] "stop_gained,inframe_insertion"                                
#> [26] "stop_gained,NMD_transcript_variant"                           
#> [27] "stop_gained,splice_region_variant"                            
#> [28] "start_lost"                                                   
#> [29] "start_lost,splice_region_variant"                             
#> [30] "start_lost,start_retained_variant"                            
#> [31] "stop_lost"                                                    
#> [32] "stop_lost,splice_region_variant"                              
#> [33] "stop_retained_variant"                                        
```
