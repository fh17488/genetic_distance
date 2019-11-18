%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Commands for Task 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fetching the mitochondrial genome (mtDNA) of the animal 'Ovis aries' commonly known as sheep from GenBank.
sheepmtdna = getgenbank('AF010406','SequenceOnly',true);
%Fetching the number of elements and the number of bytes:
whos sheepmtdna
%Fetching base count
bases=basecount(sheepmtdna)
%Fetching the local base composition plot as well as the CG-content plot using a sliding window of size = 830
%The content plot with a window size of 830 was not reported.
%window size 1: 830 = (num elements 16616 / 20).
%bases=ntdensity(sheepmtdna,'window',830)
%Fetching the local base composition plot as well as the CG-content plot using a sliding window of size = 183. The following content plot has been reported.
%window size 2: 183 = (num elements 16616 / 91).
bases=ntdensity(sheepmtdna,'window',183)
%For a correct window size a CG-plot will show an anomoly in the CG frequency iff there has been a horizontal gene transfer into the genome of the
%species. This means that a genetic sequence was inserted into the genome of the species from another organisim that does not belong to the taxonomical hierarchy of the given species, in this case the sheep.
%However, only basic statistics that describe the data have been asked and therefore the aforementioned analysis not requried and hence will not be done.
% NOT REQUIRED... next step:
%The dimer frequency and motif bias is shown below. Matlab command: [dimers,matrix] = dimercount(sheepmtdna)
%[dimers,matrix] = dimercount(sheepmtdna)
%Signficiant ORFs with atleast 46 codons is:
orf_46=seqshoworfs(sheepmtdna,'GeneticCode',1,'MinimumLength',46);
%Fetching the mitochondrial genome (mtDNA) of the animal 'Ovis Aries' commonly known as sheep from GenBank.
OAries_mtDNA = getgenbank('AF010406','SequenceOnly',true);
%Translation of Cytochrome B of 'Ovis Aries':
OAries_CYTB_ORF_start_loc = 14159;
OAries_CYTB_ORF_stop_loc = 15297;
OAries_CYTB_AA = nt2aa(OAries_mtDNA(OAries_CYTB_ORF_start_loc:OAries_CYTB_ORF_stop_loc),'GeneticCode',2,'ACGTOnly','FALSE')
%To validate download the CYTB_AA sequence from GenBank using following command:
%OAries_CYTB_AA_DB = getgenpept('AKP54989','SequenceOnly',true)
%Translation of Cytochrome C sub-unit 1 of 'Ovis Aries':
OAries_CYTC_SU1_ORF_start_loc = 5332;
OAries_CYTC_SU1_ORF_stop_loc = 6875;
OAries_CYTC_SU1_AA = nt2aa(OAries_mtDNA(OAries_CYTC_SU1_ORF_start_loc:OAries_CYTC_SU1_ORF_stop_loc),'GeneticCode',2,'ACGTOnly','FALSE')
%Translation of Cytochrome C sub-unit 2 of 'Ovis Aries':
OAries_CYTC_SU2_ORF_start_loc = 7019;
OAries_CYTC_SU2_ORF_stop_loc = 7701;
OAries_CYTC_SU2_AA = nt2aa(OAries_mtDNA(OAries_CYTC_SU2_ORF_start_loc:OAries_CYTC_SU2_ORF_stop_loc),'GeneticCode',2,'ACGTOnly','FALSE')
%Translation of Cytochrome C sub-unit 3 of 'Ovis Aries':
OAries_CYTC_SU3_ORF_start_loc = 8616;
OAries_CYTC_SU3_ORF_stop_loc = 9398;
OAries_CYTC_SU3_AA = nt2aa(OAries_mtDNA(OAries_CYTC_SU3_ORF_start_loc:OAries_CYTC_SU3_ORF_stop_loc),'GeneticCode',2,'ACGTOnly','FALSE')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 3: Bulding the Phylogenetic Tree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Retrieving DNA sequences for the CYTB coding region in the mtDNA of each of the 5 species:
%Cytochrome B of 'Ovis Aries':
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE = getgenbank('KP229289','SequenceOnly',true);
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_ARIES_CYTB_DNA_SEQ = OVIS_ARIES_CYTB_DNA_SEQ_SOURCE(OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Cytochrome B of 'Ovis Vignei Blanfordi':
OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE = getgenbank('EU366050','SequenceOnly',true);
OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ = OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE(OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Translation of Cytochrome B of 'Ovis Ammon Hodgsoni':
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE = getgenbank('JX673911','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Translation of Cytochrome B of 'Ovis Nivicola':
OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE = getgenbank('AJ867265','SequenceOnly',true);
OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_NIVICOLA_CYTB_DNA_SEQ = OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE(OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_NIVICOLA_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Translation of Cytochrome B of 'Ovis Orientalis Gmelini':
OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE = getgenbank('FJ936203','SequenceOnly',true);
OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ = OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE(OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Translation of Cytochrome B of 'Ovis Canadensis':
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14156
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15294
OVIS_CANADENSIS_CYTB_DNA_SEQ = OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE = getgenbank('DQ093614','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
CAPRA_HIRCUS_CYTB_DNA_SEQ = CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Creating structure with sequence data of original species Ovis Aries and 5 nearby species found on BLAST.
%In order to make the tree a rooted tree the Goat species has been added as an outgroup:
%        Species Description    GenBank Accession
data = {'Sheep'					OVIS_ARIES_CYTB_DNA_SEQ;
        'Blanford Urial'		OVIS_VIGNEI_BLANFORDI_CYTB_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ;
        'Snow Sheep'			OVIS_NIVICOLA_CYTB_DNA_SEQ;
        'Armenian Mouflon'		OVIS_ORIENTALIS_GMELINI_CYTB_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_DNA_SEQ;
		'Goat'					CAPRA_HIRCUS_CYTB_DNA_SEQ;
       };

for ind = 1:length(data)
    related_species_of_sheep(ind).Header   = data{ind,1};
    related_species_of_sheep(ind).Sequence = data{ind,2};
end

distances = seqpdist(related_species_of_sheep,'Method','Jukes-Cantor','Alpha','DNA');
NJtree = seqneighjoin(distances,'equivar',related_species_of_sheep)
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree of Related Species of Sheep using Jukes-Cantor model');
ylabel('Evolutionary distance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 4.1: To calculate the pair-wise distance between species based on their CYTB protein-coding DNA sequences. 
%The Jukes-Cantor model has been used to correct the estimate of the genetic distance.
%The result is generated as a matrix. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE = getgenbank('KP229289','SequenceOnly',true);
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_ARIES_CYTB_DNA_SEQ = OVIS_ARIES_CYTB_DNA_SEQ_SOURCE(OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE = getgenbank('JX673911','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14156
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15294
OVIS_CANADENSIS_CYTB_DNA_SEQ = OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE = getgenbank('DQ093614','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
CAPRA_HIRCUS_CYTB_DNA_SEQ = CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14233
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15371
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

%Creating structure with sequence data of 5 species in order to calculate the pair-wise distance between each nucleotide sequence using the Jukes-Cantor model:
data = {'Sheep'					OVIS_ARIES_CYTB_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTB_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTB_DNA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTB_DNA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

D = seqpdist(CYTB_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alphabet','NT')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 4.2: To calculate the pair-wise distance between species based on their CYTC subunit 1 protein-coding DNA sequences. 
%The Jukes-Cantor model has been used to correct the estimate of the genetic distance.
%The result is generated as a matrix. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KP981380','SequenceOnly',true);
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX - OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ = OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JX101654','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5333
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6876
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX - OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5330
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6873
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX - OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ = OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('GU229280','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX - CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ = CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5406
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6949
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

%Creating structure with sequence data of 5 species in order to calculate the pair-wise distance between each nucleotide sequence using the Jukes-Cantor model:
data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

D = seqpdist(CYTC_UNIT1_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alphabet','NT')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 5.1: Generating a Phylogenetic Tree inclusive of all species
%Data: CYTB DNA sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE = getgenbank('KP229289','SequenceOnly',true);
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_ARIES_CYTB_DNA_SEQ = OVIS_ARIES_CYTB_DNA_SEQ_SOURCE(OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE = getgenbank('JX673911','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14156
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15294
OVIS_CANADENSIS_CYTB_DNA_SEQ = OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE = getgenbank('DQ093614','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
CAPRA_HIRCUS_CYTB_DNA_SEQ = CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14233
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15371
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE = getgenbank('AF084078','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ = SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE(SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

data = {'Sheep'					OVIS_ARIES_CYTB_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTB_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTB_DNA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTB_DNA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

distances = seqpdist(CYTB_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alpha','DNA');
NJtree = seqneighjoin(distances,'equivar',CYTB_DNA_SEQS_VECTOR)
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree with outgroup - CYTB DNA sequences - using Jukes-Cantor model');
ylabel('Evolutionary distance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 5.2: Generating a Phylogenetic Tree inclusive of all species
%Data: CYTC subunit 1 DNA sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KP981380','SequenceOnly',true);
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ = OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JX101654','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5333
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6876
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5330
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6873
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ = OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('GU229280','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ = CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5406
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6949
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JF681040','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5353
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6911
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ = SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE(SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

distances = seqpdist(CYTC_UNIT1_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alpha','DNA');
NJtree = seqneighjoin(distances,'equivar',CYTC_UNIT1_DNA_SEQS_VECTOR)
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree with outgroup - CYTC Subunit 1 DNA sequences - using Jukes-Cantor model');
ylabel('Evolutionary distance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 5.3: Generating a Phylogenetic Tree inclusive of all species
%Data: CYTB Amino Acid Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTB_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTB_AA_SEQ = getgenpept('AAD54455','SequenceOnly',true);

data = {'Sheep'					OVIS_ARIES_CYTB_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_AA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTB_AA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTB_AA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTB_AA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTB_AA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

distances = seqpdist(CYTB_AA_SEQS_VECTOR,'Method','Jukes-Cantor','Alphabet','AA');
NJtree = seqneighjoin(distances,'equivar',CYTB_AA_SEQS_VECTOR)
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree with outgroup - CYTB Amino Acid sequences - using Jukes-Cantor model');
ylabel('Evolutionary distance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 5.4: Generating a Phylogenetic Tree inclusive of all species
%Data: CYTC subunit 1 Amino Acid Sequences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTC_UNIT1_AA_SEQ = getgenpept('AKT25940','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ = getgenpept('AFN06203','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ = getgenpept('AEK79957','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ = getgenpept('ACZ61158','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ = getgenpept('ANF99669','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTC_UNIT1_AA_SEQ = getgenpept('AEG23805','SequenceOnly',true);

data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTC_UNIT1_AA_SEQ;
       };
	   
for ind = 1:length(data)
    CYTC_UNIT1_AA_SEQS_VECTOR(ind).Header   = data{ind,1};
    CYTC_UNIT1_AA_SEQS_VECTOR(ind).Sequence = data{ind,2};
end

distances = seqpdist(CYTC_UNIT1_AA_SEQS_VECTOR,'Method','Jukes-Cantor','Alphabet','AA');
NJtree = seqneighjoin(distances,'equivar',CYTC_UNIT1_AA_SEQS_VECTOR)
h = plot(NJtree,'orient','top');
title('Neighbor-Joining Distance Tree with outgroup - CYTC subunit 1 Amino Acid sequences - using Jukes-Cantor');
ylabel('Evolutionary distance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Task 5.5: Generating a Consensus Tree
%Data: Weighted average of genetic distance from DNA coding sequences of CYTB and CYTC subunit 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE = getgenbank('KP229289','SequenceOnly',true);
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_ARIES_CYTB_DNA_SEQ = OVIS_ARIES_CYTB_DNA_SEQ_SOURCE(OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE = getgenbank('JX673911','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14156
OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15294
OVIS_CANADENSIS_CYTB_DNA_SEQ = OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE = getgenbank('DQ093614','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
CAPRA_HIRCUS_CYTB_DNA_SEQ = CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 14233
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 15371
CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE = getgenbank('AF084078','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_START_INDEX = 1
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX = 1139
SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ = SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE(SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_START_INDEX:SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ_SOURCE_STOP_INDEX);

CYTB_DNA_SEQ_DATA = {'Sheep'					OVIS_ARIES_CYTB_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTB_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTB_DNA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTB_DNA_SEQ;
       };
	   
for ind = 1:length(CYTB_DNA_SEQ_DATA)
    CYTB_DNA_SEQS_VECTOR(ind).Header   = CYTB_DNA_SEQ_DATA{ind,1};
    CYTB_DNA_SEQS_VECTOR(ind).Sequence = CYTB_DNA_SEQ_DATA{ind,2};
end

DISTANCE_FROM_CYTB_DNA_SEQS = seqpdist(CYTB_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alpha','DNA');

OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KP981380','SequenceOnly',true);
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
OVIS_ARIES_CYTC_UNIT1_DNA_SEQ = OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_ARIES_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JX101654','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5333
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6876
OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ = OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JN181255','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5330
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6873
OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ = OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE(OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('GU229280','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5329
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6872
CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ = CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('KU942399','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5406
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6949
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ = CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE(CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE = getgenbank('JF681040','SequenceOnly',true);
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX = 5353
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX = 6911
SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ = SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE(SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_START_INDEX:SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ_SOURCE_STOP_INDEX);

CYTC_UNIT1_DNA_SEQ_DATA = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_DNA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_DNA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_DNA_SEQ;
        'Goat'					CAPRA_HIRCUS_CYTC_UNIT1_DNA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_DNA_SEQ;
		'Gray Dolphin'			SOTALIA_FLUVIATILIS_CYTC_UNIT1_DNA_SEQ;
       };
	   
for ind = 1:length(CYTC_UNIT1_DNA_SEQ_DATA)
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Header   = CYTC_UNIT1_DNA_SEQ_DATA{ind,1};
    CYTC_UNIT1_DNA_SEQS_VECTOR(ind).Sequence = CYTC_UNIT1_DNA_SEQ_DATA{ind,2};
end

DISTANCE_FROM_CYTC_UNIT1_DNA_SEQS = seqpdist(CYTC_UNIT1_DNA_SEQS_VECTOR,'Method','Jukes-Cantor','Alpha','DNA');

LABELS = {'Sheep'					
        'Tibetan Argali'		
        'Bighorn Sheep'			
        'Goat'					
        'Alashan Red Deer'		
		'Gray Dolphin'			
       };

for ind = 1:length(LABELS)
    LABELS_VECTOR(ind).Header   = LABELS{ind,1};
end


weights = [sum(DISTANCE_FROM_CYTB_DNA_SEQS) sum(DISTANCE_FROM_CYTC_UNIT1_DNA_SEQS)];
weights = weights / sum(weights);
WEIGHTED_Avg_Distance = DISTANCE_FROM_CYTB_DNA_SEQS .* weights(1) + DISTANCE_FROM_CYTC_UNIT1_DNA_SEQS .* weights(2);
CONSENSUS_TREE = seqlinkage(WEIGHTED_Avg_Distance,'average',LABELS_VECTOR);
plot(CONSENSUS_TREE,'type','angular');
title('Consensus Tree from CYTB and CYTC Subunit 1 DNA sequences');
%Task 7.1: Perform multiple alignment and find polymorphic sites
%Data = CYTB Amino Acid Sequences
%Sheep and related species without Goat and without Deer.
OVIS_ARIES_CYTB_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);

sheep_neighbours_only_data = {'Sheep'					OVIS_ARIES_CYTB_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_AA_SEQ;        
       };
	   
for ind = 1:length(sheep_neighbours_only_data)
    CYTB_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY(ind).Header   = sheep_neighbours_only_data{ind,1};
    CYTB_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY(ind).Sequence = sheep_neighbours_only_data{ind,2};
end

seqalignviewer(multialign(CYTB_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY),'Alphabet','AA')
%Task 7.2: Perform multiple alignment and find polymorphic sites
%Data = CYTB Amino Acid Sequences
%Sheep and related species with Goat and without Deer.
OVIS_ARIES_CYTB_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);

Goat_only_data = {'Sheep'					OVIS_ARIES_CYTB_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_AA_SEQ;    
		'Goat'					CAPRA_HIRCUS_CYTB_AA_SEQ;
       };
	   
for ind = 1:length(Goat_only_data)
    CYTB_AA_SEQS_VECTOR_WITH_GOAT(ind).Header   = Goat_only_data{ind,1};
    CYTB_AA_SEQS_VECTOR_WITH_GOAT(ind).Sequence = Goat_only_data{ind,2};
end

seqalignviewer(multialign(CYTB_AA_SEQS_VECTOR_WITH_GOAT),'Alphabet','AA')
%Task 7.3: Perform multiple alignment and find polymorphic sites
%Data = CYTB Amino Acid Sequences
%Sheep and related species with Deer and without Goat.
OVIS_ARIES_CYTB_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTB_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTB_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTB_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);

Deer_only_data = {'Sheep'					OVIS_ARIES_CYTB_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTB_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTB_AA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTB_AA_SEQ;
       };
	   
for ind = 1:length(Deer_only_data)
    CYTB_AA_SEQS_VECTOR_WITH_DEER(ind).Header   = Deer_only_data{ind,1};
    CYTB_AA_SEQS_VECTOR_WITH_DEER(ind).Sequence = Deer_only_data{ind,2};
end

seqalignviewer(multialign(CYTB_AA_SEQS_VECTOR_WITH_DEER),'Alphabet','AA')
%Task 7.4: Perform multiple alignment and find polymorphic sites
%Data = CYTC UNIT 1 Amino Acid Sequences
%Sheep and related species without Goat and without Deer.
OVIS_ARIES_CYTC_UNIT1_AA_SEQ = getgenpept('AKT25940','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ = getgenpept('AFN06203','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ = getgenpept('AEK79957','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ = getgenpept('ACZ61158','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ = getgenpept('ANF99669','SequenceOnly',true);

sheep_neighbours_only_data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ;        
       };
	   
for ind = 1:length(sheep_neighbours_only_data)
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY(ind).Header   = sheep_neighbours_only_data{ind,1};
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY(ind).Sequence = sheep_neighbours_only_data{ind,2};
end

seqalignviewer(multialign(CYTC_UNIT1_AA_SEQS_VECTOR_WITH_NEIGHBOURS_ONLY),'Alphabet','AA')
%Task 7.5: Perform multiple alignment and find polymorphic sites
%Data = CYTC UNIT 1 Amino Acid Sequences
%Sheep and related species with Goat and without Deer.
OVIS_ARIES_CYTC_UNIT1_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);

Goat_only_data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ;    
		'Goat'					CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ;
       };
	   
for ind = 1:length(Goat_only_data)
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_GOAT(ind).Header   = Goat_only_data{ind,1};
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_GOAT(ind).Sequence = Goat_only_data{ind,2};
end

seqalignviewer(multialign(CYTC_UNIT1_AA_SEQS_VECTOR_WITH_GOAT),'Alphabet','AA')
%Task 7.6: Perform multiple alignment and find polymorphic sites
%Data = CYTC UNIT 1 Amino Acid Sequences
%Sheep and related species with Deer and without Goat.
OVIS_ARIES_CYTC_UNIT1_AA_SEQ = getgenpept('AKP54989','SequenceOnly',true);
OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ = getgenpept('AGK62343','SequenceOnly',true);
OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ = getgenpept('AEK79967','SequenceOnly',true);
CAPRA_HIRCUS_CYTC_UNIT1_AA_SEQ = getgenpept('AAY98750','SequenceOnly',true);
CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ = getgenpept('ANF99679','SequenceOnly',true);

Deer_only_data = {'Sheep'					OVIS_ARIES_CYTC_UNIT1_AA_SEQ;
        'Tibetan Argali'		OVIS_AMMON_HODGSONI_CYTC_UNIT1_AA_SEQ;
        'Bighorn Sheep'			OVIS_CANADENSIS_CYTC_UNIT1_AA_SEQ;
        'Alashan Red Deer'		CERVUS_ELAPHUS_ALXAICUS_CYTC_UNIT1_AA_SEQ;
       };
	   
for ind = 1:length(Deer_only_data)
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_DEER(ind).Header   = Deer_only_data{ind,1};
    CYTC_UNIT1_AA_SEQS_VECTOR_WITH_DEER(ind).Sequence = Deer_only_data{ind,2};
end

seqalignviewer(multialign(CYTC_UNIT1_AA_SEQS_VECTOR_WITH_DEER),'Alphabet','AA')
