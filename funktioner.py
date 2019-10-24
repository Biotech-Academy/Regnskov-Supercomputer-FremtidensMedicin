#!/usr/bin/env python3

#import nødvendige pakker
from IPython.display import Image
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import requests
import json
from bioservices import KEGG
k = KEGG(verbose=False)
from dna_features_viewer import BiopythonTranslator
from dna_features_viewer import GraphicFeature, GraphicRecord, CircularGraphicRecord   
from IPython.display import YouTubeVideo
from IPython.display import HTML
from google.colab import files
from IPython.display import Image

#Definere funktioner
def Hent_information_fra_GenBank_file(information):
    fil_navn = '/root/Colab_Notebooks_BA_Bioinformatik_Regnskov_Projekt/Exercises/Artemisia annua.gb'
    record = SeqIO.read(fil_navn, "genbank")
    #record.annotations.keys()

    #Hvilken slags data er det?
    if information == "data_type":
        print(record.annotations['molecule_type'])

    #Hvilken type sekventerings type er dataen lavet med?
    elif information == "sekventerings_type":
        print(record.annotations['keywords'])

    #Hvad er taxonomyen?
    elif information == "taxonomy":
        print(record.annotations['taxonomy'])

    #Print al information
    elif information == "al_information":
        print(record)

    #Print DNA 
    elif information == "dna":
        print(record.seq[0:10000])
        
    else:
        print("Du har ikke indtastet informations type korrekt")
        
    return  
  
class ChangeFeatures(BiopythonTranslator):
    def compute_feature_color(self, feature):
        if feature.type == "CDS":
            return "blue"
        elif feature.type == "assembly_gap":
            return "green"
        else:
            return "white"

    def compute_feature_label(self, feature):
        if feature.type == 'assembly_gap':
            return "Gap"
        elif feature.type == "CDS":
            return "Gen"
        #else:
        #    return BiopythonTranslator.compute_feature_label(feature)

    def compute_filtered_features(self, features):
        #"""Vis ikke andet end CDS og gaps"""
        qualifier = "locus_tag" #codon_start, locus_tag, protein_id
        #feature_type = "CDS" #mRNA, gene, gap_assembly...

        return [
            feature for feature in features
            if (feature.type != "CDS") or (feature.type != "assembly_gap")
            or ("CT" in str(feature.qualifiers[qualifier]))]
        
def Visualiser_sekvens(fil):
    graphic_record = ChangeFeatures().translate_record(fil)
    ax, _ = graphic_record.plot(figure_width=20)
    ax.figure.tight_layout()
    return

def display_yotube_video(url, **kwargs):
    id_ = url.split("=")[-1]
    return YouTubeVideo(id_, **kwargs)
  
def Hent_genetisk_informatik_API_funktion(art,protein,antal):    
    #response = requests.get("http://string-db.org/api/json/resolve?identifier=PfEMP1&species=Plasmodium falciparum")
    string = "http://string-db.org/api/json/resolve?identifier="+protein+"&species="+art
    response = requests.get(string)
    #print(response.status_code)
    #print(response.content)

    #Gem i json format
    data = response.json()
    print("Antal hits:", len(data), "\n")
    
    if antal != "alt":
        print("Se første", antal, "hits:", "\n")
        [print(i) for i in data[0:antal]]
        return
        
    #print("\n","Se alle IDs:", "\n")
    #print([i["stringId"] for i in data])
    else:
        #data_2 = [i["stringId"] for i in data]
        return data #_2  
      
def Visualisér_sekvens(Gener):
    Gener = 0
    features=[
        GraphicFeature(start=6258, end=7320, strand=-1, color="#ffd700",label="psbD"),
        GraphicFeature(start=7488, end=8568, strand=-1, color="#ffcccc",label="psbA"),
        GraphicFeature(start=81944, end=82268, strand=-1, color="#cffccc",label="PetE")
    ]
    
    record = GraphicRecord(sequence_length=92268, features=features)
    record.plot()

    features2=[
        GraphicFeature(start=6258, end=7320, strand=-1, color="#ffd700",label="psbD"),
        GraphicFeature(start=7488, end=8568, strand=-1, color="#ffcccc",label="psbA")
    ]
    record2 = GraphicRecord(sequence_length=12000, features=features2)
    record2.plot()
    
    return 


def Find_genet_funktion(Gen_navn, file_navn):
    file_navn = '/root/Colab_Notebooks_BA_Bioinformatik_Regnskov_Projekt/Exercises/Artemisia annua.gb'
  
    #Sæt faste parameter
    qualifier = "gene"
    feature_type = "CDS"
    x = 0
    
    #Find genet og print den rigtige information
    for gb_record in SeqIO.parse(open(file_navn,"r"), "genbank"):  
        for (index, feature) in enumerate(gb_record.features) :
            if feature.type==feature_type:
                if qualifier in feature.qualifiers :
                    for value in feature.qualifiers[qualifier]:
                        if value == Gen_navn:
                            print("Fundet genet:",value,2*"\n","Information om genet:","\n",feature,"\n","Index nr:",index)
                            x = 1
    if x != 1:
        print("Genet er ikke fundet")
    return 


def sekvens_transformer_funktion(fil_navn, ønsket_sekvens, index):
    for gb_record in SeqIO.parse(open(fil_navn,"r"), "genbank"):
        gb_feature = gb_record.features[index]
        coding_dna = gb_feature.extract(gb_record.seq) #Hent codende dna for gen
        template_dna = coding_dna.reverse_complement() #Se den modsatte streng
        messenger_rna = coding_dna.transcribe() #Konverter til mRNA
        protein_sequence = messenger_rna.translate() #Konverter til protein sekvens
        
        #Hm er det virkelig den bedste måde at gøre det på C, hvad med de skriver det selv?? 
        if ønsket_sekvens == "dna":
            print(coding_dna)

        elif  ønsket_sekvens == "reverse_dna":
            print(template_dna)

        elif ønsket_sekvens == "mRNA":
            print(messenger_rna)

        elif ønsket_sekvens == "protein":
            print(protein_sequence)  
            
        else:
            print("Du har ikke indtastet typen af den ønsket sekvens korrekt. Prøv venligst igen")
    return 
  
