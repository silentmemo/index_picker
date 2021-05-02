import sys 
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mode

### using the index below for development, will use dual combinatorial index

def load_all_index(filepath):
    print("Loading available index form : ",filepath)
    ### build a dictionary of available index by pandas df, then to dict
    index_dict = pd.read_csv(filepath,sep="\t").to_dict(orient= 'index')
    return index_dict

def load_chosen_index(filepath):
    print("Loading chosen index form : ",filepath)
    ### build a dictionary of available index by pandas df, then to dict
    index_dict = pd.read_csv(filepath,sep="\t").to_dict(orient= 'index')
    return index_dict

def get_sequence_from_dict_to_list(inputDict,field_name):
    outputList = []
    for item in inputDict:
        # print(inputDict[item].get("Weight"))
        outputList.append((inputDict[item].get(field_name),float(inputDict[item].get("Weight"))))
    return outputList
    # return outputList
def get_sequence_length(inputSeq):
    return len(inputSeq)
def get_weight_sum(inputWeightList):
    counter = 0
    for w in inputWeightList:
        counter = counter + w
    return counter
def get_base_composition(inputList):
    ### expecting a list of tuples of valid sequence
    #### take the first item of the input list to tell the length of the sequence
    seq_length = get_sequence_length(inputList[0][0])
    weight_list = []
    for idx,item in enumerate(inputList):
        weight = item[1]
        weight_list.append(float(weight))
    weight_sum = get_weight_sum(weight_list)
    # print(seq_length)
    #### goal is to get a per position matrix of base distribution, store that in a dictionary 
    string_composition_list =[]
    
    for base in range(seq_length):
        position_composition_dict = {"A":0, "T":0,"C": 0, "G":0}
        for idx,item in enumerate(inputList):
            sequence = item[0]
            weight = item[1] / weight_sum
            position_composition_dict[sequence[base]] = position_composition_dict[sequence[base]] + weight
            
        string_composition_list.append(position_composition_dict)
    return string_composition_list
def update_weight(inputDict,new_weight):
    # print(inputDict)
    remaining_weight = 1 - new_weight
    for item in inputDict:
        inputDict[item]["Weight"] = inputDict[item]["Weight"] * remaining_weight
        print(inputDict[item]["Weight"])
def get_mode_weight(inputList):
    weight_list = []
    for item in inputList:
        weight = item[1]
        weight_list.append(weight)
    return float(mode(weight_list)[0])
def get_variance(inputList):
    outputList = []
    for item in inputList:
        outputList.append(np.var(item))
    return sum(outputList)
def main():

    print("Picking indexes")
    ava_index_dict = load_all_index(sys.argv[1])
    # ava_index_dict = load_all_index("CD_index.txt")
    chosen_index_dict = load_chosen_index(sys.argv[2])
    # chosen_index_dict = load_chosen_index("chosen_indexes.txt")
    another_flag = True
    while another_flag != False:
        chosen_I7_seq_list = get_sequence_from_dict_to_list(chosen_index_dict,"I7_sequence")
        chosen_I5_seq_list = get_sequence_from_dict_to_list(chosen_index_dict,"I5_sequence")

        ### get the base composition from chosen indexes
        # print(chosen_I7_seq_list)
        I7_base_composition = get_base_composition(chosen_I7_seq_list)
        I5_base_composition = get_base_composition(chosen_I5_seq_list)
        I7_df  = pd.DataFrame(I7_base_composition)
        I5_df  = pd.DataFrame(I5_base_composition)
        # print(chosen_I7_seq_list)

        ### with the base composition found, design a maximization function such that it maximize the base balance 
        ### the base balance is maximized when the varience between base are close to 0 

        I7_values = I7_df.values.tolist()
        chosen_I7_var = get_variance(I7_values)
        
        ### brute force to select the next best I7 by looping through the list 
        #### need to assume the weight of the upcoming index, choose mode of chosen list
        mode_chosen_I7_weight = get_mode_weight(chosen_I7_seq_list)
        I7_test_result_dict = {}
        for item in ava_index_dict:
            ### copy a new temp list for each simulation
            temp_chosen = chosen_I7_seq_list.copy()
            if ava_index_dict[item]["index_type"] == "I7":
                ### let assume to put this index into chosen, then calculate the variance again 
                ava_sequence  = ava_index_dict[item]["sequence"]
                ava_index_name = ava_index_dict[item]["index_name"]
                
                ava_tuple = (ava_sequence, mode_chosen_I7_weight)
                temp_chosen.append(ava_tuple)
                ### compute varience again 
                temp_base_composition = get_base_composition(temp_chosen)
                temp_df = pd.DataFrame(temp_base_composition)
                temp_values = temp_df.values.tolist()
                temp_var = get_variance(temp_values)
                
                ### put the result in the dictionary 
                I7_test_result_dict[ava_index_name] = temp_var

        #### do it for I5 below
        I5_values = I5_df.values.tolist()
        chosen_I5_var = get_variance(I5_values)
        

        mode_chosen_I5_weight = get_mode_weight(chosen_I5_seq_list)
        I5_test_result_dict = {}
        for item in ava_index_dict:
            temp_chosen = chosen_I5_seq_list.copy()
            if ava_index_dict[item]["index_type"] == "I5":
                ava_sequence  = ava_index_dict[item]["sequence"]
                ava_index_name = ava_index_dict[item]["index_name"]
                ava_tuple = (ava_sequence, mode_chosen_I5_weight)
                temp_chosen.append(ava_tuple)
                temp_base_composition = get_base_composition(temp_chosen)
                temp_df = pd.DataFrame(temp_base_composition)
                temp_values = temp_df.values.tolist()
                temp_var = get_variance(temp_values)
                I5_test_result_dict[ava_index_name] = temp_var

        #### make a combination list of I7 and I5, then rank the combined minimum 
        combo_var_list = []
        for i7 in I7_test_result_dict:
            for i5 in I5_test_result_dict:
                # print(i7,i5)
                sum_var = I7_test_result_dict[i7]+I5_test_result_dict[i5]
                combo_dict = {"I7_index_name":i7,"I5_index_name":i5,"sum_of_variance":sum_var}
                combo_var_list.append(combo_dict)
        combo_var_df = pd.DataFrame(combo_var_list)
        # print(combo_var_df.sort_values(by=['sum_of_variance']))

        ### the first minimum is the best option, if it has not been in the chosen in the chosen list
        #### select the best index
        best_index = combo_var_df.sort_values(by=['sum_of_variance'])[:1]
        ### check if it is in the chosen index
        chosen_I7_I5_list = []
        for chosen_index in chosen_index_dict:
            I7_name = chosen_index_dict[chosen_index]["I7_index_name"]
            I5_name = chosen_index_dict[chosen_index]["I5_index_name"]
            I7name_I5name = I7_name+"|"+I5_name
            chosen_I7_I5_list.append(I7name_I5name)
        best_index_I7_I5 = best_index["I7_index_name"] +"|"+ best_index["I5_index_name"]
        
        best_index_I7_I5 = best_index_I7_I5.to_string(index = False)

        ### validate that it is not a member in the chosen_I7_I5_list

        if best_index_I7_I5 in chosen_I7_I5_list:
            print(best_index_I7_I5, "Index already exist")
            ###  TODO try the next best combination
            
        else:
            print(best_index_I7_I5,"is the best option.")
            ### generate a new member in the dictionary , look up the sequence from ava
            ava_df = pd.DataFrame.from_dict(ava_index_dict,orient="index")
            new_I7_name = best_index_I7_I5.split("|")[0].strip()
            new_I5_name = best_index_I7_I5.split("|")[1].strip()
            new_I7_sequence = ava_df[ava_df["index_name"] == new_I7_name]["sequence"].to_string(index = False).strip()
            new_I5_sequence = ava_df[ava_df["index_name"] == new_I5_name]["sequence"].to_string(index = False).strip()
            
            user_input_weight = input("Enter weight (fraction) of new library: ")
            user_input_weight = float(user_input_weight)
            update_weight(chosen_index_dict,user_input_weight)
            new_index_dict_value = {'I7_index_type': 'I7', 'I7_index_name': new_I7_name, 'I7_sequence': new_I7_sequence, 'I5_index_type': 'I5', 'I5_index_name': new_I5_name, 'I5_sequence': new_I5_sequence, 'Weight': user_input_weight}
            new_index_dict_key = len(chosen_index_dict)
            chosen_index_dict[new_index_dict_key] = new_index_dict_value
            # print(chosen_index_dict)
            #### ask if user want to add another index 
            user_input_want_more_index = input("Want to add another index? (Y/N): ")
            
            if user_input_want_more_index == "Y":
                another_flag = True
            else:
                another_flag = False
                print (pd.DataFrame.from_dict(chosen_index_dict,orient="index"))
                i7_fig = I7_df.plot(kind="bar", stacked=True).get_figure()
                i7_fig.savefig('i7_base_composition.png')
                i5_fig = I5_df.plot(kind="bar", stacked=True).get_figure()
                i5_fig.savefig('i5_base_composition.png')

main()
