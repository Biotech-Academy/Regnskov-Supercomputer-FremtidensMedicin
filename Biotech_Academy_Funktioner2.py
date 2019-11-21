# Funktioner er lavet til Biotech Academy projeketet: Regnskov + Supercomputer = Fremtidens medicin. 
# Funktionerne og materialet præsenteret er ikke akademisk korrekt, men er blevet ændret for undervisnings formål. 
# ==============================================================================


#https://widdowquinn.github.io/2018-03-06-ibioic/02-sequence_databases/09-KEGG_programming.html

# Show images inline
from IPython.display import Image

# Standard library packages
import io
import os

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

# Import Pandas, so we can use dataframes
import pandas as pd



# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# Some code to return a Pandas dataframe, given tabular text
def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

def KEGG(input1, input2):
    # Perform the query
    result = REST.kegg_info("kegg").read()

    # Print overview
    if input1 == "info" and input2 == "alt":
        return print(result)
    
    # Get all entries in the PATHWAY database as a dataframe
    elif input1 == "pathway_overview" and input2 == "alle":
        result = REST.kegg_list("pathway").read()
        return to_df(result)

    #Print alle biosynteseveje
    elif input1 == "print_pathway":
        if input2 == "alle":
              result = REST.kegg_get("map01100", "image").read()
              img = Image(result, width=1500, height=1000)
        else:
              result = REST.kegg_get(input2, "image").read()
              img = Image(result)
        return img

    #Print hvilket som helst map 
    #elif information.startswith("map") and molekyle != "alle":
    #    result = REST.kegg_get(information, "image").read()
    #    return Image(result)
    
    #Find the compund vanillin
    elif input1 == "find_molekyle" and input2 != None: 
        result = REST.kegg_find("compound",input2).read() #cpd:C00755
        return print(result)

    elif input1 == "info_molekyle" and input2 != None: #cpd:C00755
        # Get the entry information for vanillin
        result = REST.kegg_get(input2).read()
        return print(result)

    # Display molekylær struktur for cpd:C00051 (vanillin)
    elif input1 ==  "molekyle billede" and input2 != None::
        result = REST.kegg_get(input2, "image").read() #"cpd:C00755"
        return Image(result)

    elif input1 == "Enzyme molekyle" and input2 != None::
        result = REST.kegg_find("enzyme", input2).read()
        return to_df(result)
    

