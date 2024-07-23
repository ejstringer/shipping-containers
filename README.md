# shipping-containers
Pest insects in shipping containers eDNA and eRNA anlaysis

## Goals

* Difference between eDNA and eRNA
* Detection ~ time since arrival
* Detection ~ country (number/different container histories)
* Risk assessment (visual detection and food/habitat of pest in container)


## Metadata

Original data files can be found in the folder labelled *Metabarcoding publication 2023* (owed by Alejandro Trujillo-Gonzalez)

**Shipping container meta data:**

*Metabarcoding publication 2023/Container history/*

* __Additional diagnostics extract.xlsx__ (Shipping containers)


Notes:

Loading country

*	This is derived from multiple sources including the container event data provided by the shipping companies.
*	Indicates the country in which a container was loaded onto a vessel for export to another country. The container could be full or empty. 
*	Unfortunately, there are some gaps and errors in the load country data. This is the result of stitching together event data from different systems/sources.
*	Please note, where you see duplicate loading countries for the same date and different destination countries – this is a likely indication of a **transhipment event**.

Origin country  

*	This is derived from ICS data and therefore only existing for events with a destination of Australia i.e. when the container was imported into Australia. 
*	It indicates the origin of the goods in the container i.e. most likely where they were packed into the container. 

Goods description

*	This is derived from ICS, AIMS data and therefore only existing for events with a destination of Australia i.e. when the container was imported into Australia.

Shipping company goods description 

*	Goods description provided by the shipping company. 

Destination country    

*	Generally the country where a container is discharged from a vessel.
*	Can include transshipment countries. 

Arrival sequence number and Is latest arrival

*	Arrival sequence number 1 = the last voyage for the container i.e. when it was imported to Australia prior to unpacking and sampling. 
*	Ordering of the arrival sequence number allows you to track the movement history of the container. 
*	Unfortunately, the sequencing of transshipment events appears to be jumbled in some cases. 



**Species specific files:**

*Metabarcoding publication 2023/Nucleic acid data/*

* __khapra beetle_2023.xlsx__      (Khapra beetle,*Trogoderma granarium*) 
* __Electric Ant_2023.csv__        (Electric Ant, *Wasmannia auropunctata*)
* __Spotted lantern fly_2023.csv__ (Spotted lantern fly, *Lycorma delicatula*)
* __BMSB_2023.csv__                (Brown marmorated stink bug, *Halyomorpha halys*)  
* __Asian spongy moth_2023.csv__   (Asian spongy moth, *Lymantria dispar asiatica*) 

**Metabarcoding:**

*Metabarcoding publication 2023/Nucleic acid data/Metabarcoding/*

* __ASV_all_tax_count.tsv__  (eDNA)
* __ASV_cdna_tax_count.tsv__ (eRNA)


## Code


File 1:

* File name: filter_metabarcodingDNA_file-Arthropoda.R
* Purpose: filtering the very large ASV file with all phylum to only arthropoda asvs 

File 2:

* File name: initial_data_exploration.R
* Purpose: To identify the container ids to link all the datasets by

File 3:




