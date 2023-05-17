# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:46:34 2023

@author: Sahil
"""

import streamlit as st
from selenium import webdriver
driver = webdriver.Chrome(executable_path = 'C:/Users/Sahil/.wdm/drivers/chromedriver/win32/113.0.5672/chromedriver.exe')
from selenium.webdriver.common.keys import Keys
from selenium.common.exceptions import NoSuchElementException
from time import sleep
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
import pandas as pd
import time
import regex as re
from PIL import Image
import matplotlib.pyplot as plt
import seaborn as sns

class IEDB_MHC1_Immunogenicity:
    def __init__(self):
        self.base_url = 'http://tools.iedb.org/immunogenicity/'
      
    def get_data(self, source_filepath):
        epit_list = []
        length_list = []
        immuno_score_list = []
    #home page
        driver.get(self.base_url)
        sleep(3)
    #choose file
        choose_file = driver.find_element_by_xpath('//*[@id="id_sequence_file"]')
        choose_file.send_keys(source_filepath)
    #submit button
        submit = driver.find_element_by_xpath('//*[@id="input-form"]/table/tbody/tr[6]/th/div/input[1]')
        submit.click()
    #wait until output table appears
        WebDriverWait(driver, 5).until(EC.visibility_of_element_located((By.XPATH, '//*[@id="result_table"]')))
    #finding all elements from the table
        epit = driver.find_elements_by_class_name('sequence')
        length = driver.find_elements_by_xpath('//*[@id="result_table"]/tbody/tr/td[2]')  
        score = driver.find_elements_by_xpath('//*[@id="result_table"]/tbody/tr/td[3]')
    #adding them into list
        for ep, leng, sco in zip(epit, length, score):
            epit_list.append(ep.text)
            length_list.append(int(leng.text))
            immuno_score_list.append(float(sco.text))
   #close the window     
        #driver.close()
    #return dataframe
        df = pd.DataFrame({'Epitope':epit_list, 'Length':length_list, 'Immunogenicity':immuno_score_list})
        return df
    
    def save_results(self, df, dest_filepath):
        df.to_csv(dest_filepath, index = False)
        
        
class Vaxijen_Antigenicity:
    def __init__(self):
        self.base_url = 'http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
        
    def get_data(self, source_filepath):
        epit_list = []
        score_list = []
        is_antigen = []
        #home page
        driver.get(self.base_url)
        sleep(3)
        #choosing file
        choose_file =  driver.find_element_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[1]/td[2]/p/input')
        choose_file.send_keys(source_filepath)
        #selecting the target organism
        target_org = driver.find_element_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[2]/td[2]/p/select/option[2]')
        target_org.click()
        #submit button
        submit = driver.find_element_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/form/table/tbody/tr[3]/td[2]/input[1]')
        submit.click()
        WebDriverWait(driver, 5).until(EC.visibility_of_element_located((By.XPATH, '/html/body/div/table/tbody/tr[3]/td[3]/font/b')))
        epit = driver.find_elements_by_xpath('//font[@face = "Courier New"]')
        score_tag = driver.find_elements_by_xpath('/html/body/div/table/tbody/tr[4]/td[3]/table/tbody/tr/td/b')
        anti = driver.find_elements_by_xpath('//td/b/font[@color = "#0000FF"]')
        #score is in <b> tag. There are variuos <b> tags. So we will use regex to check
        regex = '[+-]?[0-9]+\.[0-9]+'
        for i in range(0, len(score_tag)):
            if re.search(regex, score_tag[i].text):
                score_list.append(float(score_tag[i].text))  
        for epi, is_ant in zip(epit, anti):
            epit_list.append(epi.text)
            is_antigen.append(is_ant.text)
        #driver.close()
        df = pd.DataFrame({"Epitope":epit_list, "Antigenicity":score_list, "Is_Antigen":is_antigen})
        return df
    
    def save_results(self, df, dest_filepath):
        df.to_csv(dest_filepath)  
        
        
class ToxinPred_Toxicity:
    def __init__(self):
        self.base_url = 'https://webs.iiitd.edu.in/raghava/toxinpred2/batch.html'
        
    def get_data(self, source_filepath):
        driver.get(self.base_url)
        sequences = []
        toxicity_scores = []
        toxic = []
        chose_file = driver.find_element_by_xpath('/html/body/header/div[2]/section/form/table/tbody/tr/td/font/p/font[2]/input')
        chose_file.send_keys(source_filepath)
        sleep(2)
        submit = driver.find_element_by_xpath('/html/body/header/div[2]/section/form/table/tbody/tr/td/font/font/p[3]/font/font[2]/input[2]')
        submit.click()
        WebDriverWait(driver, 60).until(EC.visibility_of_element_located((By.XPATH, '/html/body/header/div[2]/main/div/table[2]')))
        seq = driver.find_elements_by_xpath('/html/body/header/div[2]/main/div/table[2]/tbody/tr/td[1]')
        score = driver.find_elements_by_xpath('/html/body/header/div[2]/main/div/table[2]/tbody/tr/td[5]')
        tox = driver.find_elements_by_xpath('/html/body/header/div[2]/main/div/table[2]/tbody/tr/td[6]')
        for s, t, to in zip(seq, score, tox):
            sequences.append(s.text)
            toxicity_scores.append(float(t.text))
            toxic.append(to.text)
        #driver.close()
        df1 = pd.DataFrame({'Sequence':sequences, 'Toxicity_score':toxicity_scores, 'Toxicity':toxic})
        return df1
    
    def convert_fasta_to_dataframe(self, source_filepath):
        with open(source_filepath, 'r') as f:
            pep = f.readlines()
        seq = []
        epit = []
        for i in range(0, len(pep)):
            if i%2 == 0:
                seq.append(pep[i].split('>')[1].split('\n')[0])
            else:
                epit.append(pep[i].split('\n')[0])   
        df2 = pd.DataFrame({'Sequence':seq, 'Epitope':epit})
        return df2
    
    def merge_results(self, df1, df2):
        df = pd.merge(df1, df2, on = 'Sequence')
        return df
    
    def save_results(self, df, dest_filepath):
        df.to_csv(dest_filepath, index = False)
        
class AllerTop_Allergenicity:
    def __init__(self):
        self.base_url = 'https://www.ddg-pharmfac.net/AllerTOP/index.html'
        
    def get_data(self, source_filepath):
        with open(source_filepath, 'r') as f:
            pep = f.readlines()
        epitope_sequences = []
        for i in range(len(pep)):
            epitope_sequences.append(pep[i].split('\n')[0])
        driver.get(self.base_url)    
        epit_list = []
        allergenicity = []
        for i in range(0, len(epitope_sequences)):
            #textbox for entering sequence
            epit_seq = driver.find_element_by_xpath('//*[@id="sequence"]')
            epit_seq.send_keys(epitope_sequences[i])
            #pressing the submit button
            submit = driver.find_element_by_xpath('//*[@id="protein_sequence"]/table/tbody/tr[3]/td[1]/input')
            submit.click()
            sleep(3)
            #peptide name
            epit_list.append(epitope_sequences[i])
            #allergenicity
            is_aller = driver.find_element_by_xpath('//*[@id="box"]/h4[2]').text
            allergenicity.append(is_aller)
            #coming back to home page
            driver.get(self.base_url)
            sleep(3)
        #close the window
        driver.close()
        dataframe = pd.DataFrame({'Epitope':epit_list, 'Allergenicity':allergenicity})
        return dataframe
    
    def save_results(self, data, dest_filepath):
        data.to_csv(dest_filepath, index = False)


st.title('Belyntic precision vaccine GmbH') 
image = Image.open(r'C:\Users\Sahil\.spyder-py3\belyntic_logo.jpg') 
st.image(image)      
st.write('Epitopes score')



source_filepath_csv = r'C:\Users\Sahil\Omdena\Development of Vaccines for JC Virus\9per_peptides.csv'
source_filepath_fasta = r'C:\Users\Sahil\Omdena\Development of Vaccines for JC Virus\mhc1_epitopes_fasta.txt'
df = pd.read_csv(source_filepath_csv)


if st.sidebar.button('Run all pipeline'):
    st.info('Currently running IEDB Immunogenicity')
    imm = IEDB_MHC1_Immunogenicity()
    i_data = imm.get_data(source_filepath_csv)
    st.write(i_data.head())
    st.info('Immunogenicity calculated successfully')

    st.info('Currently running Vaxigen Antigenicity')
    va = Vaxijen_Antigenicity()
    a_data = va.get_data(source_filepath_fasta)
    st.write(a_data.head())
    st.info('Antigenicity calculated successfully')

    st.info('Currently running Toxinpred Toxicity')
    tt = ToxinPred_Toxicity()
    data1 = tt.get_data(source_filepath_fasta)
    data2 = tt.convert_fasta_to_dataframe(source_filepath_fasta)
    t_data = tt.merge_results(data1, data2)
    st.write(t_data.head())
    st.info('Toxicity calculated successfully')
    
    st.info('Currently running Allertop Allergenicity')
    aa = AllerTop_Allergenicity()
    al_data = aa.get_data(source_filepath_csv)
    st.write(al_data.head())
    st.info('Allergenicity calculated successfully')


    st.info('Merging datasets')
    final_a = pd.merge(i_data, a_data, on = 'Epitope')
    final_b = pd.merge(final_a, t_data, on = 'Epitope')
    final = pd.merge(final_b, al_data, on = 'Epitope')
    st.write(final)

    @st.cache
    def convert_df(df):
        return df.to_csv(index = False)
    final_csv = convert_df(final)
    st.download_button('Download data as csv', data = final_csv, file_name= 'epitopes_score.csv')


indi_tool = st.sidebar.selectbox('Run individual tools', ['None',  'IEDB', 'Vaxigen', 'Toxinpred', 'Allertop'])
if indi_tool == 'IEDB':
    st.info('Currently running IEDB Immunogenicity')
    imm = IEDB_MHC1_Immunogenicity()
    i_data = imm.get_data(source_filepath_csv)
    st.write(i_data)
    st.success('Immunogenicity calculated successfully')

if indi_tool == 'Vaxigen':
    st.info('Currently running Vaxigen Antigenicity')
    va = Vaxijen_Antigenicity()
    a_data = va.get_data(source_filepath_fasta)
    st.write(a_data)
    st.success('Antigenicity calculated successfully')

if indi_tool == 'Toxinpred':
    st.info('Currently running Toxinpred Toxicity')
    tt = ToxinPred_Toxicity()
    data1 = tt.get_data(source_filepath_fasta)
    data2 = tt.convert_fasta_to_dataframe(source_filepath_fasta)
    t_data = tt.merge_results(data1, data2)
    st.write(t_data)
    st.success('Toxicity calculated successfully')
    
if indi_tool == 'Allertop':
    st.info('Currently running Allertop Allergenicity')
    aa = AllerTop_Allergenicity()
    al_data = aa.get_data(source_filepath_csv)
    st.write(al_data)
    st.success('Allergenicity calculated successfully')


file = st.sidebar.file_uploader('Upload csv file', type = 'csv')
if file is not None:
    df = pd.read_csv(file)
    st.write(df)
    Wimm = st.number_input('Enter immunogenicity weight', min_value=0.0, max_value=1.0, step = 0.01)
    Want = st.number_input('Enter antigenicity weight', min_value=0.0, max_value=1.0, step = 0.01)
    Wtox = st.number_input('Enter toxicity weight', min_value=0.0, max_value=1.0, step = 0.01)
    df['rank'] = Wimm*df['Immunogenicity'] + Want*df['Antigenicity'] - Wtox*df['Toxicity_score']
    df.sort_values(by = 'rank', ascending = False, inplace = True)
    st.write(df)
    st.info('Rank calculated successfully')
     
    if st.button("Heatmap"):
        correlation = df[['Immunogenicity', 'Antigenicity', 'Toxicity_score', 'rank']].corr()
        fig = plt.figure(figsize = (6,6))
        sns.heatmap(data = correlation, cmap = 'coolwarm', annot = True, fmt = '.1f')
        st.pyplot(fig = fig)

    if st.button('Scatterplot'):
        rank_formula_features = ['Immunogenicity', 'Antigenicity', 'Toxicity_score']
        fig, axes = plt.subplots(1,3, figsize = (20, 10)) 
        for param, ax in zip(rank_formula_features, axes.flatten()):
            sns.scatterplot(data = df, x =param, y = 'rank', ax = ax)
        st.pyplot(fig=fig, use_container_width = False)

    if st.button('Paramter contribution plot'):
        fig = plt.figure(figsize=(15, 8))
        sns.scatterplot(data=df[['Immunogenicity', 'Antigenicity', 'Toxicity_score', 'rank']], s=100)
        plt.xticks(rotation=90)
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        st.pyplot(fig=fig)




