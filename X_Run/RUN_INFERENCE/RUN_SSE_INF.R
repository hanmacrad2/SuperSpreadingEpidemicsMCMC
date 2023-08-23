#************************
#* IMPLEMENT
#***********************

#FOLDER
folder_type = 'data_sse'
CURRENT_FOLDER = GET_FOLDER_TIME_STAMP(folder_type, array_index)
create_folder(CURRENT_FOLDER)

#RUN_INFERENCE_SSE
df_sse_inf = RUN_INFERENCE_SSE(DATA_FOLDER)
