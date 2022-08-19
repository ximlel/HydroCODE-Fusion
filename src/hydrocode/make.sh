#!/bin/bash

### Compile the program
# make clean
make


### Run the program
## GRP_Book
:<<!
 ./hydrocode.out GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 1     LAG
 ./hydrocode.out GRP_Book/6_1_LAG GRP_Book/6_1_LAG 1 2_GRP LAG
 ./hydrocode.out GRP_Book/6_1_EUL GRP_Book/6_1_EUL 1 1     EUL
 ./hydrocode.out GRP_Book/6_1_EUL GRP_Book/6_1_EUL 1 2_GRP EUL
 ./hydrocode.out GRP_Book/6_2_1   GRP_Book/6_2_1   1 1     EUL
 ./hydrocode.out GRP_Book/6_2_1   GRP_Book/6_2_1   1 2_GRP EUL
 ./hydrocode.out GRP_Book/6_2_3   GRP_Book/6_2_3   1 1     EUL
 ./hydrocode.out GRP_Book/6_2_3   GRP_Book/6_2_3   1 2_GRP EUL
!

## GRP_direct
:<<!
 ./hydrocode.out GRP_direct/9_1_a GRP_direct/9_1_a 1 1     EUL
 ./hydrocode.out GRP_direct/9_1_a GRP_direct/9_1_a 1 2_GRP EUL 
 ./hydrocode.out GRP_direct/9_1_b GRP_direct/9_1_b 1 1     EUL
 ./hydrocode.out GRP_direct/9_1_b GRP_direct/9_1_b 1 2_GRP EUL
 ./hydrocode.out GRP_direct/9_1_d GRP_direct/9_1_d 1 1     EUL
 ./hydrocode.out GRP_direct/9_1_d GRP_direct/9_1_d 1 2_GRP EUL
 ./hydrocode.out GRP_direct/9_1_e GRP_direct/9_1_e 1 1     EUL
 ./hydrocode.out GRP_direct/9_1_e GRP_direct/9_1_e 1 2_GRP EUL
!

## Artifi_Heat_Conduct
:<<!
 ./hydrocode.out Artifi_Heat_Conduct/4_22 Artifi_Heat_Conduct/4_22 1 1     LAG
 ./hydrocode.out Artifi_Heat_Conduct/4_22 Artifi_Heat_Conduct/4_22 1 2_GRP LAG 
 ./hydrocode.out Artifi_Heat_Conduct/4_2 Artifi_Heat_Conduct/4_2 1 1     LAG
 ./hydrocode.out Artifi_Heat_Conduct/4_2 Artifi_Heat_Conduct/4_2 1 2_GRP LAG 
 ./hydrocode.out Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 1 1     LAG
 ./hydrocode.out Artifi_Heat_Conduct/4_3 Artifi_Heat_Conduct/4_3 1 2_GRP LAG 
 ./hydrocode.out Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 1 1     LAG
 ./hydrocode.out Artifi_Heat_Conduct/4_5 Artifi_Heat_Conduct/4_5 1 2_GRP LAG
!
