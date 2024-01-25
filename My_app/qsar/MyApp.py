import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import numpy as np

# Molecular descriptor calculator


def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download


def filedownload(df):
    csv = df.to_csv(index=False)
    # strings <-> bytes conversions
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building


def build_model(input_data, smiles, verbose=False):
    st.header("Calcule Lipinski Discriptors and predict the bioactivity")
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    # st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='Active_Inactive')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df_ppred = pd.concat([molecule_name, prediction_output], axis=1)
    ################################

    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if (i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, descriptors], axis=1)
    df_d = pd.DataFrame()
    df_drug = pd.DataFrame(["Respect"])
    df_ndrug = pd.DataFrame(["dont Respect"])

    for i in range(df.shape[0]):
        if df.MW[i] < 500 and df.LogP[i] <= 5 and df.NumHDonors[i] <= 5 and df.NumHAcceptors[i] <= 10:
            df_d = df_d._append(df_drug, ignore_index=True)
        else:
            df_d = df_d._append(df_ndrug, ignore_index=True)

    df_d = df_d.rename(columns={0: "lipinski_rules"})
    df_f = pd.concat([df, df_d], axis=1)
    df_final = pd.concat(
        [df_f, df_ppred['Active_Inactive']], axis=1)
    st.write(df_final)
    st.markdown(filedownload(df_final), unsafe_allow_html=True)

    # return df
    # st.write(df)
    # st.markdown(filedownload(df), unsafe_allow_html=True)


# calcule Lipinski Disc

# def lipinski_calc(smiles, verbose=False):

    # return df_f


###########################################################

    # Logo image
image = Image.open("logo.png")

st.sidebar.image(image, use_column_width=False)

imagee = Image.open("logo2.png")

st.image(imagee, use_column_width=True, caption='acetylcholinesterase protein')

# Page title
st.markdown("""
# Predicting Acetylcholinesterase Inhibition and Lipinski Rule Evaluation


---
""")

# Sidebar
with st.header('Upload your txt data'):
    uploaded_file = st.file_uploader(
        "File with two columns: molecule_id, molecule_smile", type=['txt'])


if st.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi', sep='\t', header=False, index=False)

    st.header('**Original input data**')
    st.write(load_data)

    with st.spinner("Calculating descriptors..."):
        desc_calc()

    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('descriptors_output.csv')
    st.write(desc)
    st.write(desc.shape)

    # Read descriptor list used in previously built model
    # st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    # st.write(desc_subset)
    # st.write(desc_subset.shape)

    # Apply trained model to make prediction on query compounds
    with st.spinner("Lipinski and classification..."):
        df_pred = build_model(desc_subset, load_data[0])

    # Lipinski disc

    # lipinski = lipinski_calc(load_data[0])
    # df_final = pd.concat(
    #    [lipinski, df_pred['class|Active|-|Inactive|']], axis=1)
    # st.write(df_final)
    # st.markdown(filedownload(df_final), unsafe_allow_html=True)
else:
    st.info('Upload input data and click predict.')
