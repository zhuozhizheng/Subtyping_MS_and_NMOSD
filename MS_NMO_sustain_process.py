from ntpath import join
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sqlalchemy import true
import pySuStaIn
import statsmodels.formula.api as smf
from scipy import stats
from scipy.stats import ttest_ind
import sklearn.model_selection
df1=pd.read_csv('data/PFDR_value.csv')
#df=pd.read_csv('data/SuStain_MS_type.csv')
#df=pd.read_csv('data/tem_clinical_data_Volume_20221119_complete_SUStainHC.csv')
#df=pd.read_csv('data/R_ALL_Sustain1.csv')
#df=pd.read_csv('data/lv_wmh/Data_list_Cluster_number5_new.csv')
df=pd.read_csv('data/Sustain_adjusted_data.csv')
#data_type='AQP4_Pos'
data_type='MS'
#data_type='MOGAD'
#data_type='AQP4Pos_NMOSD'
#data_type='AQP4Neg_NMOSD'
#data_type='AQP4Neg_NMOSD+AQP4Pos_NMOSD'
#data_type='AQP4Neg_NMOSD+AQP4Pos_NMOSD+MOGAD+MS'
#data_type='all'
#data_type='AQP4_Neg'
#data_type='AQP4_Neg+AQP4_Pos'
N_S_max = 5 
N_startpoints = 10 
output_folder = os.path.join(os.getcwd(),'data', 'WorkshopOutput_1212_nmo')
def ttest_arr(train_x,train_y):
    train_0=[]
    train_1=[]
    for i in range(len(train_x)):
        if train_x[i]=='':
            continue
        if train_y[i]<1:
            train_0.append(train_x[i])
        else:
            train_1.append(train_x[i])
    return train_0,train_1
def get_top10(df,class_index,biomarkers):
    p_dict={}
    p_value=[]
    y=[0 if t=='HC' else 1 for t in df[class_index]]
    for t in biomarkers:
        t_0,t_1=ttest_arr(df[t],y)
        t_r=ttest_ind(t_0,t_1)
        p_value.append(t_r[1])
        p_dict[t]=t_r[1]
    p_value=sorted(p_value)
    if len(p_value)>10:
        out_biomarkers=[t for t in biomarkers if p_dict[t]<p_value[10]]
    else:
        out_biomarkers=[biomarkers]
    print(out_biomarkers)
    print(p_value)
    print(p_dict)
    return(out_biomarkers)

def get_biomakers(df=df1):
    biomakers = [ df.values[i,0] for i in range(len(df)) if df[data_type][i]<0.05 ]
    return biomakers

def eval(data=df,biomarkers=None,class_index='Group'):
    zdata=pd.read_csv('data/Sustain_adjusted_data_follow_up.csv')
    print(len(zdata))
    if biomarkers==None:
        biomarkers =[t for t in data.columns[-19:]]
    # for each biomarker cal z-score
    norm_index=data[class_index]=='HC'
    for biomarker in biomarkers:
        #print(biomarker)
        t_data=data[biomarker][norm_index]
        norm_mean=np.mean(t_data)
        norm_std=np.std(data[biomarker][norm_index])
        z_score = -(zdata.loc[:,biomarker] - norm_mean) / norm_std
        zdata.loc[:,biomarker] = z_score
    N = len(biomarkers)         # number of biomarkers
    #print(len(zdata))

    #new_dict={}
    #ms_index=np.zeros_like(zdata[class_index])
    #ms_index=ms_index>0
    #for k in zdata.keys():
    #    new_dict[k]=zdata[k][ms_index]
    #zdata=pd.DataFrame(new_dict)
    return zdata

def test(data=df,biomarkers=None,class_index='Group'):
    #data.Diagnosis.value_counts()
    if biomarkers==None:
        #biomarkers =get_top10(data,data.columns[0],data.columns[1:])
        biomarkers =[t for t in data.columns[-19:]]
        #biomarkers = get_biomakers()
    print(biomarkers)
    if class_index==None:
        class_index=data.columns[-1]
    print(class_index)
    zdata = pd.DataFrame(data,copy=True)

    # for each biomarker cal z-score
    norm_index=zdata[class_index]=='HC'
    for biomarker in biomarkers:
        t_data=zdata[biomarker][norm_index]
        norm_mean=np.mean(t_data)
        norm_std=np.std(zdata[biomarker][norm_index])
        z_score = -(data.loc[:,biomarker] - norm_mean) / norm_std
        zdata.loc[:,biomarker] = z_score
    N = len(biomarkers)         # number of biomarkers
    new_dict={}
    #ms_index=zdata[class_index]=='AQP4_Pos'
    ms_index=np.zeros_like(zdata[class_index])
    if data_type=='all':
        ms_index=zdata[class_index]!='HC'
    for d in data_type.split('+'):
        ms_index+=zdata[class_index]==d
    ms_index=ms_index>0
    for k in zdata.keys():
        new_dict[k]=zdata[k][ms_index]
    zdata=pd.DataFrame(new_dict)

    SuStaInLabels = biomarkers
    Z_vals = np.array([[1,2,3]]*N)     # Z-scores for each biomarker
    Z_max  = np.array([5]*N)           # maximum z-score
    # Input the settings for z-score SuStaIn
    # To make the tutorial run faster I've set 
    # N_startpoints = 10 and N_iterations_MCMC = int(1e4)
    # I recommend using N_startpoints = 25 and 
    # N_iterations_MCMC = int(1e5) or int(1e6) in general though

    N_iterations_MCMC = int(1e4)
    #dataset_name = 'AQP4_Pos'
    dataset_name = data_type
    print(len(zdata))

    # Initiate the SuStaIn object
    sustain_input = pySuStaIn.ZscoreSustain(
                                  zdata[biomarkers].values,
                                  Z_vals,
                                  Z_max,
                                  SuStaInLabels,
                                  N_startpoints,
                                  N_S_max, 
                                  N_iterations_MCMC, 
                                  output_folder, 
                                  dataset_name, 
                                  True)


    # # Run SuStaIn!
    # make the output directory if it's not already created
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)


    zdata.to_csv(os.path.join(output_folder,'ori'+data_type+'test.csv'),index=False)
    samples_sequence,   \
    samples_f,          \
    ml_subtype,         \
    prob_ml_subtype,    \
    ml_stage,           \
    prob_ml_stage,      \
    prob_subtype_stage  = sustain_input.run_sustain_algorithm()
    print(len(data))
    e_data=eval(data)
    print(len(e_data))
    e_data.to_csv(os.path.join(output_folder,'ori'+data_type+'eval.csv'),index=False)
    ml_subtype, prob_ml_subtype, ml_stage, prob_ml_stage, prob_subtype, prob_stage, prob_subtype_stage=sustain_input.subtype_and_stage_individuals_newData(e_data[biomarkers].values,samples_sequence,samples_f,len(e_data))
    e_data['ml_subtype']=ml_subtype
    e_data['prob_ml_subtype']=prob_ml_subtype
    e_data['ml_stage']=ml_stage
    e_data['prob_ml_stage']=prob_ml_stage
    e_data.to_csv(os.path.join(output_folder,'new'+data_type+'eval.csv'),index=False)


def cross_validation():
    cvics=[]
    for i in range(10):
        t=cross_validation_fold(i)
        print(t)
        cvics.append(t)
    np.savetxt(os.path.join(output_folder,'cvic_10.txt'),np.array(cvics))
def cross_validation_fold(fold):
    class_index='Group'
    zdata=pd.read_csv(os.path.join(output_folder,'ori'+data_type+'test.csv'))
    biomarkers =[t for t in zdata.columns[-19:]]
    N = len(biomarkers)         # number of biomarkers
    
    SuStaInLabels = biomarkers
    Z_vals = np.array([[1,2,3]]*N)     # Z-scores for each biomarker
    Z_max  = np.array([5]*N)           # maximum z-score
    # Input the settings for z-score SuStaIn
    # To make the tutorial run faster I've set 
    # N_startpoints = 10 and N_iterations_MCMC = int(1e4)
    # I recommend using N_startpoints = 25 and 
    # N_iterations_MCMC = int(1e5) or int(1e6) in general though

    N_iterations_MCMC = int(1e4)
    #dataset_name = 'AQP4_Pos'
    dataset_name = data_type
    if not os.path.isdir(os.path.join(output_folder,str(fold))):
        os.mkdir(os.path.join(output_folder,str(fold)))

    # Initiate the SuStaIn object
    sustain_input = pySuStaIn.ZscoreSustain(
                                  zdata[biomarkers].values,
                                  Z_vals,
                                  Z_max,
                                  SuStaInLabels,
                                  N_startpoints,
                                  N_S_max, 
                                  N_iterations_MCMC, 
                                  os.path.join(output_folder,str(fold)),
                                  dataset_name, 
                                  True)


    # # Run SuStaIn!
    # make the output directory if it's not already created
    # # Cross-validation

    # Finally, it is often difficult to decide how many subtypes best fit your data. This question should ideally be evaluated using cross-validation. This way, the likelihood metrics are generated for data that the model has not yet seen. 
    # 
    # FYI we may not have nough time to run this full cross-validation during the workshop, but it's good for you to run it yourself. SuStaIn also support parallelized cross-validation!

    # choose the number of folds - here i've used three for speed but i recommend 10 typically
    
    N_folds = 10

    # generate stratified cross-validation training and test set splits
    labels = zdata[class_index]
    cv = sklearn.model_selection.StratifiedKFold(n_splits=N_folds, shuffle=True)
    cv_it = cv.split(zdata, labels)

    # SuStaIn currently accepts ragged arrays, which will raise problems in the future.
    # We'll have to update this in the future, but this will have to do for now
    test_idxs = []
    for train, test in cv_it:
        test_idxs.append(test)
    test_idxs = np.array(test_idxs,dtype='object')


    # In[ ]:


    # perform cross-validation and output the cross-validation information criterion and
    # log-likelihood on the test set for each subtypes model and fold combination
    CVIC, loglike_matrix     = sustain_input.cross_validate_sustain_model(test_idxs)
    return CVIC


    # # Choosing the optimal number of subtypes
    # The optimal number of subtypes is chosen using the CVIC, shown below. The CVIC is an information criterion (like the AIC/BIC/WAIC) that balances model complexity with model accuracy, with a lower CVIC indicating a better balance between the two. Generally speaking, the model with the lowest CVIC is the best. However, you do sometimes get a very small improvement (less than ~6) in the CVIC with a more complex model, in which case I would tend to favour the less complex (i.e. fewer subtypes) model.
    # 
    # Another useful metric to look at is the log-likelihood of each subtypes model on the test set, also shown below. A better model should show a consistent improvement in the test set log-likelihood across folds.

    # In[ ]:


    # go through each subtypes model and plot the log-likelihood on the test set and the CVIC
    #print("CVIC for each subtype model: " + str(CVIC))
    #print("Average test set log-likelihood for each subtype model: " + str(np.mean(loglike_matrix, 0)))

    #plt.figure(0)    
    #plt.plot(np.arange(N_S_max,dtype=int),CVIC)
    #plt.xticks(np.arange(N_S_max,dtype=int))
    #plt.ylabel('CVIC')  
    #plt.xlabel('Subtypes model') 
    #plt.title('CVIC')

    #plt.figure(1)
    ##df_loglike = pd.DataFrame(data = loglike_matrix, columns = ["s_" + str(i) for i in range(sustain_input.N_S_max)])
    ##df_loglike.boxplot(grid=False)
    #plt.ylabel('Log likelihood')  
    #plt.xlabel('Subtypes model') 
    #plt.title('Test set log-likelihood across folds')


    ## Another useful output of the cross-validation that you can look at are positional variance diagrams averaged across cross-validation folds. These give you an idea of the variability in the progression patterns across different training datasets.


    ##this part estimates cross-validated positional variance diagrams
    #for i in range(N_S_max):
    #    sustain_input.combine_cross_validated_sequences(i+1, N_folds)



    ##N_S_selected = 2

    ##pySuStaIn.ZscoreSustain._plot_sustain_model(sustain_input,samples_sequence,samples_f,M,subtype_order=(0,1))
    ##_ = plt.suptitle('SuStaIn output')

    ##sustain_input.combine_cross_validated_sequences(N_S_selected, N_folds)
    ##_ = plt.suptitle('Cross-validated SuStaIn output')

    #plt.show()





if __name__=='__main__':
    test()
    cross_validation()

