import requests
import hashlib
import shutil
import os
import subprocess
import datetime
import time
import shlex
import csv
import json
from pprint import pprint

now = datetime.datetime.now

#Check if Windows or Linux
#print(os.name)
if os.name == 'nt': #Windows
	gdc_location = "C:\\Users\\Deitrickc\\Downloads\\Genomic Programs\\gdc-client"
	gdc_program  = os.path.join(gdc_location, 'gdc-client.exe')
	histology_filename = "C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Clinical Data\\histology_diagnoses.txt"
else:
	histology_filename = "/home/upmc/Documents/Variant_Discovery_Pipeline/api_files/histology_diagnoses.txt"
	gdc_location = "/home/upmc/Programs/gdc_data_transfer_tool"
	local_file_api_filename = "/home/upmc/Documents/Variant_Discovery_Pipeline/api_files/local_file_api.json" #Path to a local copy of the file api
	local_case_api_filename = "/home/upmc/Documents/Variant_Discovery_Pipeline/api_files/local_case_api.json" #Path to a local copy of the case api
	gdc_program  = os.path.join(gdc_location, 'gdc-client-2016-10-19')
user_token = os.path.join(gdc_location, 'tokens', max(list(os.listdir(os.path.join(gdc_location, "tokens")))))

print(user_token)



def timeout_command(command, timeout):
    """call shell-command and either return its output or kill it
    if it doesn't normally exit within timeout seconds and return None"""
    import subprocess, datetime, os, time, signal
    start = datetime.datetime.now()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
      time.sleep(0.1)
      now = datetime.datetime.now()
      if (now - start).seconds> timeout:
        os.kill(process.pid, signal.SIGKILL)
        os.waitpid(-1, os.WNOHANG)
        return None
    return process.stdout.read()

def Terminal(command, label = None, show_output = False, timeout = None):
	""" Calls the system shell """
	terminal_log = "terminal_log.txt"
	label = "{0} {1}".format(label if label is not None else "", now().isoformat())
	terminal_label  = '--'*15 + label + '--'*15 + '\n'
	terminal_label += '..'*40 + '\n'
	terminal_label += command + '\n'
	terminal_label += '..'*40 + '\n'
	terminal_label += '--'*40 + '\n'

	#Try using exceptions to catch timeout errors

	if show_output and False:
		process = os.system(command)
	else:
		#from subprocess import STDOUT, check_output
		#output = subprocess.check_output(command, stderr=subprocess.STDOUT, timeout=timeout)
		command = shlex.split(command)
		print(command)
		
		#process = subprocess.call(command, timeout=7200, shell=False)
		#process = subprocess.run(command, timeout = 7200, shell = False)
		process = subprocess.Popen(command, shell = False)
		print("Communicating with subprocess...")
		try:
			process.communicate(timeout = 3600+1800)
			process.kill()
		except:
			process.kill()

	return process

def _find_best_pair(files):
    """ If a case has more than one WXS normal/tumor read available,
        this selects the best one.
        Parameters
        ----------
            files: list<dict>
                A list containing the metadata of the files to parse.
    """
    if isinstance(files[0], str):
        files = [case_api(f, expand = 'Full') for f in files]
    normals = list()
    tumors = list()
    
    for f in files:
        if f['file_info']['tissue_type'] == 'normal':
            normals.append(f)
        elif f['file_info']['tissue_type'] == 'tumor':
            tumors.append(f)
            
    normal = sorted(normals, key = lambda s: s['file_info']['sample_id'])[0]
    tumor = sorted(tumors,   key = lambda s: s['file_info']['sample_id'])[0]

    return {'normal': normal, 'tumor': tumor}

def _get_exome_targets(catalog_number):
	""" Maps an exome target url to a local file """
	if catalog_number is not None:
		exome_targets_folder = "/home/upmc/Documents/Reference/Exome_Targets/"
		_931070 = os.path.join(exome_targets_folder , "whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed")
		_6465692001 = os.path.join(exome_targets_folder , "SeqCap_EZ_Exome_v3_capture.bed")
		catalog_nnumber = '|'.join([i for i in catalog_number.split('|') if i != ['NA']])

		local_files = {
			#Custom V2 Exome Bait, 48 RXN X 16 tubes
			'931070': _931070,
			#Nimblegen SeqCap EZ Human Exome Library v3.0
			'06465692001': _6465692001
		}
		catalog_number = max(catalog_number.split('|'), key = lambda s: len(s))
		targets = local_files.get(catalog_number)
	else:
		targets = None

	return targets

def request_file(file_id, folder = None, tries = 10):
	""" Downloads a file to the current working directory.
		Parameters
		----------
			file_id: string
				The file id of the file to download.
		Returns
		----------
			downloaded_directory: string [PATH]
				The path to the directory holding the downloaded files.
				Usually organized as 
				CWD
					File_id
						logs
							log1
							log2
						file1
						file2
		Options
		-------

		"""
	file_info = file_api(file_id)
	file_name = file_info['basic_info']['file_name']
	if folder:
		output_folder = os.path.join(folder, file_id)
	else:
		output_folder = os.path.join(os.getcwd(), file_id)
	main_output_file = os.path.join(output_folder, file_name)

	gdc_command = "{gdc} download {file_id} --token-file {token} --no-annotations" #+ "--verbose"
	if folder:
		gdc_command = gdc_command + " --dir {folder}"
		
	gdc_command = gdc_command.format(
		gdc = gdc_program,
		file_id = file_id,
		token = user_token,
		folder = folder)

	attempts = 0
	while attempts < tries:
		attempts += 1
		print("\t\t\t Attempt {0} of {1}: ".format(attempts, tries), end = ' ', flush = True)
		result = Terminal(gdc_command, show_output = False, label = 'Request File', timeout = 3600 + 1800)

		status = _verify_file(main_output_file, file_info)
		if status: break
		time.sleep(10)

	response = {
		'status': status,
		'output folder': output_folder,
		'path': main_output_file
	}

	return response

def _request(endpoint, uuid, parameters):
	url = 'https://gdc-api.nci.nih.gov/{endpoint}/{uuid}'.format(
		endpoint = endpoint,
		uuid = uuid)
	response = requests.get(url, parameters)
	#print(response.url)
	response = response.json()
	
	if 'data' in response:
		response = response['data']
	return response

def _verify_file(path, file_info):
	file_exists = os.path.isfile(path)
	if file_exists:
		file_md5sum = hashlib.md5(open(path,'rb').read()).hexdigest()
		md5_check = file_md5sum == file_info['md5sum']
	else:
		md5_check = False

	return file_exists and md5_check

def download_file(file_id, destination, basename = None):
	""" Downloads the file with the 'file_id' to the location referenced by 'path'
		Parameters
		----------
			file_id: string
				The file id of the file to download.
			destination: string
				The folder to move the downloaded files to.
				Only the files will be moved, not the folder they were
				originally downloaded in.
			basename: string
				This will replace the original filenames of the files.
				The original file extensions will be kept.
		Returns
		----------
			output_folder: string [PATH]
				The folder the files were moved to.

	"""
	cwd = os.getcwd()
	#------------------- Download the Files -------------------------
	response = request_file(file_id, destination, tries = 10)
	download_directory = response['output folder']

	#-------------- Move the files to the destination ---------------
	#Changed
	#The files are now downloaded to the correct folder via the gdc-client.
	
	if response['status'] and basename:
		for fn in os.listdir(download_directory):
			current_path = os.path.join(download_directory, fn)

			if os.path.isfile(current_path):
				#Generate the new path to the file
				current_basename, current_ext = os.path.splitext(fn)
				if basename is not None:
					current_destination = os.path.join(destination, basename + current_ext)
				else:
					current_destination = os.path.join(destination, fn)
				#Move the file, creating any required directories on the way
				if not os.path.isdir(destination):
					os.makedirs(destination)
				shutil.move(current_path, current_destination)
				
			else:
				#The path points to a directory
				if basename is not None:
					current_destination = os.path.join(destination, basename)
				else:
					current_destination = os.path.join(destination)
		#Remove the now-empty downloaded folder
		shutil.rmtree(download_directory)
		response['output folder'] = os.path.dirname(current_destination)
		response['path'] = current_destination
	else:
		pass

	return response

def _get_histology():
	"""Attempts to extract histology data from the original spreadsheets
	"""

	#histology_filename = "C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Clinical Data\\histology_diagnoses.txt"
	try:
		with open(histology_filename, 'r') as file1:
			reader = [l for l in csv.reader(file1, delimiter = '\t')]
		reader = [(i[0].lower(), i[1]) for i in reader]
		reader = dict(reader)
	except:
		harcoded_histology = """0500F1A6-A528-43F3-B035-12D3B7C99C0F	Esophagus Adenocarcinoma, NOS
								70084008-697D-442D-8F74-C12F8F598570	Esophagus Adenocarcinoma, NOS
								606DC5B8-7625-42A6-A936-504EF25623A4	Esophagus Adenocarcinoma, NOS
								CEAF98F8-517E-457A-BF29-ACFE22893D49	Esophagus Adenocarcinoma, NOS
								EE47CD59-C8D8-4B1E-96DB-91C679E4106F	Esophagus Adenocarcinoma, NOS
								61DF8A4B-95F8-40AB-A252-D00C4300C290	Esophagus Adenocarcinoma, NOS
								AB2755D2-5BD9-4E2F-B255-E813AFC8D268	Esophagus Adenocarcinoma, NOS
								6584FA8F-B4F4-4844-A71E-3BFB731AD445	Esophagus Adenocarcinoma, NOS
								078DC5F1-D2DC-4408-B785-27703C7813F1	Esophagus Adenocarcinoma, NOS
								CA400EA1-3E4E-437E-A54C-446431B741DA	Esophagus Adenocarcinoma, NOS
								F23794D7-35CC-4550-A03F-8E3CDE1A2BCD	Esophagus Adenocarcinoma, NOS
								E9AA4430-0203-413C-A8FE-BB89E95CF6A3	Esophagus Adenocarcinoma, NOS
								64AA95B2-CCF7-4BFA-8D4F-98DE09F03E71	Esophagus Adenocarcinoma, NOS
								9DFDDC49-1E66-4E26-9D71-40D748CCA0A6	Esophagus Squamous Cell Carcinoma
								7BF3F18B-A013-48D2-A743-26A57FDBE4C4	Esophagus Squamous Cell Carcinoma
								0AC08ECF-4D70-48A4-9DFD-F44032551CB9	Esophagus Squamous Cell Carcinoma
								4B8809A6-4AB3-4567-8623-37CE3FAA978C	Esophagus Squamous Cell Carcinoma
								7AB6DE57-D659-47F1-A70F-1B67650E99D7	Esophagus Squamous Cell Carcinoma
								A1A233EA-5DDA-4EEB-BD92-6667E8B3EF57	Esophagus Squamous Cell Carcinoma
								6A05288E-D051-40A5-9F4E-864AB806A5DC	Esophagus Squamous Cell Carcinoma
								6594ACD5-4D2A-4CC5-8CF8-4064E252C5F6	Esophagus Squamous Cell Carcinoma
								990C6412-6DC0-4DBA-A7C1-5BCB3751B1CC	Esophagus Adenocarcinoma, NOS
								3F7B600B-282E-4F35-9FD8-567B3FCA8273	Esophagus Squamous Cell Carcinoma
								9E54A27F-D0FA-4A5B-91BC-13E27FED2A3D	Esophagus Squamous Cell Carcinoma
								F8DBAB24-B9F4-4B8A-BFEA-57856CCF6364	Esophagus Squamous Cell Carcinoma
								A4D36E09-5D8E-4006-855C-E466682B6813	Esophagus Squamous Cell Carcinoma
								DF0125B2-5278-4DED-B3B7-9FC829FE9F63	Esophagus Squamous Cell Carcinoma
								5DA26BA5-297A-4D04-883D-32719DB3886A	Esophagus Squamous Cell Carcinoma
								41B6C417-1EC1-4DB9-91A8-11AA4DBAD05D	Esophagus Squamous Cell Carcinoma
								41B3DD5A-6413-40F5-9721-78CB218567B7	Esophagus Adenocarcinoma, NOS
								A248BFE7-59AF-43EC-8129-692F44B2B218	Esophagus Squamous Cell Carcinoma
								5E491131-B567-4E96-A2D3-69CE16CAEC50	Esophagus Squamous Cell Carcinoma
								7AA71C68-4205-4745-838B-0D00BD49A501	Esophagus Squamous Cell Carcinoma
								8F641AC4-8C31-464F-8BD1-E0661408EECC	Esophagus Adenocarcinoma, NOS
								3633AEAB-A763-438F-AD03-1B4AC8B3A621	Esophagus Squamous Cell Carcinoma
								C17331FB-C47B-40B5-8D39-9803EE37A42A	Esophagus Adenocarcinoma, NOS
								0396743B-98A6-4AD5-AA13-1B9B23E5B9F9	Esophagus Squamous Cell Carcinoma
								26E7F083-A91B-4C04-84FE-D43BE9722B18	Esophagus Squamous Cell Carcinoma
								517738F6-BA1C-43B0-8D68-81A7C34118EA	Esophagus Squamous Cell Carcinoma
								BD1C3D3B-D174-4186-BB35-364E167A1D18	Esophagus Adenocarcinoma, NOS
								97394C78-2815-4348-B4C9-DF0B3B3D6B5E	Esophagus Adenocarcinoma, NOS
								0EE03772-F38C-43C1-8484-3C4CB8459468	Esophagus Adenocarcinoma, NOS
								BCAB5381-C037-4723-904E-D01EA386B61E	Esophagus Adenocarcinoma, NOS
								7CD5ACE2-D340-487E-BEB2-1B576F6D751D	Esophagus Adenocarcinoma, NOS
								50F15B12-4DC4-403F-94B0-9D121306EE30	Esophagus Adenocarcinoma, NOS
								A6BBFC3C-BD42-47F8-AF01-99FA821EA176	Esophagus Squamous Cell Carcinoma
								EF573329-9659-4106-B391-9A33E1C732C4	Esophagus Squamous Cell Carcinoma
								A935DAEB-1EFD-4113-A614-A5E4D446C72B	Esophagus Adenocarcinoma, NOS
								6969FE5A-5993-48E5-95C5-C5C7D3D08205	Esophagus Adenocarcinoma, NOS
								06CA7FC4-4AE0-47D5-9E83-F92BA42A957B	Esophagus Squamous Cell Carcinoma
								22EE75E0-5D61-41F4-AACD-B4E09C70D2C3	Esophagus Adenocarcinoma, NOS
								209E2E16-CB88-4A8F-9CDB-10E152CAA26C	Esophagus Squamous Cell Carcinoma
								29DA08F4-C9E2-4E7A-9733-B4EC5B9C5F13	Esophagus Adenocarcinoma, NOS
								0A466142-C513-4257-85D4-4BD7CFD0EF29	Esophagus Adenocarcinoma, NOS
								0C78B1B5-F480-45AD-A253-19297BC886B4	Esophagus Adenocarcinoma, NOS
								EB4AD376-4CC1-4D2C-9B93-AD2B082C77C4	Esophagus Adenocarcinoma, NOS
								3FA7FC16-49DA-4BF6-90C4-74DF5CF68B65	Esophagus Adenocarcinoma, NOS
								CB2AF13B-CBED-4657-8FD4-FED8DB375CBA	Esophagus Adenocarcinoma, NOS
								CB419344-6148-43D4-B549-570BD09B3E05	Esophagus Adenocarcinoma, NOS
								1A578AE7-5C0B-48EF-8BAF-7189DF703A93	Esophagus Squamous Cell Carcinoma
								CABA34F2-3877-48CA-9FBC-665C3960F574	Esophagus Adenocarcinoma, NOS
								53990B20-16BD-4792-A4CA-4FD411C841F4	Esophagus Adenocarcinoma, NOS
								97031483-1185-41D6-9722-0A03B63CAA19	Esophagus Adenocarcinoma, NOS
								5530590B-8B5E-43A4-B2CA-877BFC6107F3	Esophagus Adenocarcinoma, NOS
								09EACBFD-3F4A-4D2B-B7A7-C021535E2981	Esophagus Adenocarcinoma, NOS
								D5AC647F-0C6D-40D7-87A8-FAA1CC968C82	Esophagus Adenocarcinoma, NOS
								9559F579-5D9A-4C96-9836-8A039D394795	Esophagus Adenocarcinoma, NOS
								AE6C307A-CA04-4618-B270-E8641AFD1DAA	Esophagus Adenocarcinoma, NOS
								10724958-D007-44FC-ADA2-D5874614DDBA	Esophagus Adenocarcinoma, NOS
								80BFE92B-FB00-48BB-B659-8F8EF103BC29	Esophagus Adenocarcinoma, NOS
								25EB6A8A-5BD9-480D-BAAB-6113E1B3BE51	Esophagus Squamous Cell Carcinoma
								6826BCDD-24A4-4DD1-9FAB-201CF774E166	Esophagus Adenocarcinoma, NOS
								C61536A4-761C-4CA2-A1CD-AAD738DCD940	Esophagus Adenocarcinoma, NOS
								BF3074E6-8780-490E-A8C3-4237C19CA877	Esophagus Squamous Cell Carcinoma
								FD43A93A-A26D-427B-BE04-466630F3D949	Esophagus Adenocarcinoma, NOS
								4D4E15D1-5D54-4EE3-808F-0B28DFC77871	Esophagus Squamous Cell Carcinoma
								F28F37A0-8E24-4A36-90F0-83CF5A27BE10	Esophagus Adenocarcinoma, NOS
								FB777A96-2E24-4B60-9C7E-9D5DC3E7F35A	Esophagus Adenocarcinoma, NOS
								67FAE19C-3CE9-4F91-B01E-4B130ABF632E	Esophagus Adenocarcinoma, NOS
								DF78950C-198D-4AFC-A81E-8CFA4DEA0ABD	Esophagus Adenocarcinoma, NOS
								81BE56DA-BEA2-426B-9576-950A909CC153	Esophagus Adenocarcinoma, NOS
								5AF0E222-3DC8-400B-BA61-2225921F2FD3	Esophagus Adenocarcinoma, NOS
								4B0274A5-56BD-4E69-9616-F878392EDB13	Esophagus Adenocarcinoma, NOS
								9E02A368-8858-4EC5-8C27-9E0876D387A4	Esophagus Adenocarcinoma, NOS
								5111B51D-0C7A-4B3A-847E-200B0B3EE2C4	Esophagus Squamous Cell Carcinoma
								9C32A8F1-49A4-4662-A5FB-B99F955C3B1A	Esophagus Adenocarcinoma, NOS
								A9904356-F050-4BD8-B48A-50B64F6DC3A1	Esophagus Adenocarcinoma, NOS
								90940FF1-C270-44F5-9949-F726838AAA57	Esophagus Adenocarcinoma, NOS
								19F61A54-2DB5-4E2D-8572-F294D7A41486	Esophagus Squamous Cell Carcinoma
								DA33C637-E9DB-4A67-9198-2EDA750DEAD7	Esophagus Adenocarcinoma, NOS
								318E257E-7BFF-4CA0-8F0E-ACB02516437E	Esophagus Adenocarcinoma, NOS
								66118386-6FFE-4375-97D3-01908814AE32	Esophagus Adenocarcinoma, NOS
								CF7B2DE3-F4BE-4D29-B537-499DD948C1B7	Esophagus Adenocarcinoma, NOS
								C6AE1D70-A6A2-4A46-B67A-D1BFF73A1709	Esophagus Adenocarcinoma, NOS
								D30AB67A-A40C-49B4-8FBC-D7DD5442277D	Esophagus Adenocarcinoma, NOS
								70C5CE85-4022-46C0-8BD0-F5E9009DDF41	Esophagus Squamous Cell Carcinoma
								DC4062D7-1C81-4B19-AABD-8307A0DF5029	Esophagus Adenocarcinoma, NOS
								0F06153A-0F28-4A7E-9E81-1E0B3F4412AB	Esophagus Squamous Cell Carcinoma
								AD29744C-A09D-463A-84B5-63CFEB7FF786	Esophagus Squamous Cell Carcinoma
								B3E248B5-5D40-495C-BD10-90ECF0ED6A95	Esophagus Squamous Cell Carcinoma
								AE404398-6544-4F3A-8830-44A4A227F8D5	Esophagus Squamous Cell Carcinoma
								3DF1B9DB-D900-47D2-839D-3900829CB2A9	Esophagus Squamous Cell Carcinoma
								B64A5F5A-5913-432E-9D08-0FF67A8C6B64	Esophagus Squamous Cell Carcinoma
								4C9FB311-0860-4A04-A2F6-79558AC1D98F	Esophagus Squamous Cell Carcinoma
								E9F82800-647A-4D61-BD1F-2C17C4CCAFBA	Esophagus Squamous Cell Carcinoma
								931A50FF-3149-4C4A-9456-45ED3A9DD457	Esophagus Squamous Cell Carcinoma
								0AB29DE2-EB9E-4A25-99DE-5F543CFBDB53	Esophagus Squamous Cell Carcinoma
								DAF211F2-941E-4244-8CF5-D31492269443	Esophagus Squamous Cell Carcinoma
								B05AEFC4-2D61-44B6-9451-E660D10D8A79	Esophagus Squamous Cell Carcinoma
								87FABA0E-28E7-4A7E-A7E1-12582EC68FAC	Esophagus Squamous Cell Carcinoma
								3C1F8B73-BAD7-4649-B08E-C7AFC647EDF3	Esophagus Squamous Cell Carcinoma
								A91C3F09-5F48-4B54-B669-C768F9FE9682	Esophagus Squamous Cell Carcinoma
								3C0BB5E2-AD33-4D91-A31E-2A18F1C94D2D	Esophagus Squamous Cell Carcinoma
								E89D7900-915D-4936-B1B2-CAC3670EA2F4	Esophagus Squamous Cell Carcinoma
								5080E9D1-D1E0-4CE7-A4BA-9935DD362263	Esophagus Squamous Cell Carcinoma
								FFCFA005-A04F-458E-9D1D-86143DD823E5	Esophagus Squamous Cell Carcinoma
								AA364B98-0097-43A4-976D-742FAA0276DF	Esophagus Squamous Cell Carcinoma
								C965A778-8696-47DB-A7C7-8B9F75417CBD	Esophagus Squamous Cell Carcinoma
								33C0199A-BE0C-46D1-BB1E-889BFAEDF20B	Esophagus Squamous Cell Carcinoma
								CC14603B-E32C-43F3-8218-B6ED7EE2EC18	Esophagus Squamous Cell Carcinoma
								656E3C14-1FC6-473A-813E-719753E7B408	Esophagus Squamous Cell Carcinoma
								4E7F0FFF-7B17-4DA9-95EA-F68797C427F5	Esophagus Squamous Cell Carcinoma
								D67330EE-1791-4F5A-8996-8687A2D5CCCB	Esophagus Squamous Cell Carcinoma
								2C792B77-1BF1-42F7-BD94-73C6474C6C7D	Esophagus Squamous Cell Carcinoma
								95858EA9-F71B-48D4-AA15-7DFDF21CA9E7	Esophagus Squamous Cell Carcinoma
								8C4F5492-3567-4E12-8607-2CD65FAB9BAA	Esophagus Squamous Cell Carcinoma
								7A32847F-69BB-4592-8FEB-CFD5D355ADA9	Esophagus Squamous Cell Carcinoma
								EF8BEC6F-F594-4076-9F38-157165633EF1	Esophagus Squamous Cell Carcinoma
								B782C157-6B33-4247-83DF-A82DC5E5D432	Esophagus Squamous Cell Carcinoma
								59DA377A-12A6-454A-BE01-886EDA894179	Esophagus Squamous Cell Carcinoma
								83762308-4DBA-4870-B90B-64331EBC9BAA	Esophagus Squamous Cell Carcinoma
								CAD8BD3A-FFD8-4BAF-9DE6-3AC3A1EB19E1	Esophagus Squamous Cell Carcinoma
								981066EC-ECA6-4621-8979-E935279AA1CA	Esophagus Squamous Cell Carcinoma
								37D13493-975C-432C-BD21-65F383FC66C9	Esophagus Squamous Cell Carcinoma
								CFDC8F4D-A25B-4C63-9832-58444581A505	Esophagus Squamous Cell Carcinoma
								124C1D54-1836-4F8F-920F-14047376120F	Esophagus Adenocarcinoma, NOS
								07B7DC7B-9407-4F11-9CD9-3B6DF8D1CF6C	Esophagus Squamous Cell Carcinoma
								C89F6EA5-786D-4D3C-AA93-681664A908F9	Esophagus Adenocarcinoma, NOS
								AC437816-E057-4F2C-9443-6289199431ED	Esophagus Adenocarcinoma, NOS
								DDE17514-BD86-41E7-AE6D-88DDA6D854BB	Esophagus Adenocarcinoma, NOS
								8D5788CA-9131-4604-870F-C7EEA84C9FF6	Esophagus Adenocarcinoma, NOS
								D125A297-88BE-4500-8EF8-B8BB02362F2B	Esophagus Adenocarcinoma, NOS
								7E53564C-4132-40BD-A4D5-1D35473F224F	Esophagus Adenocarcinoma, NOS
								3DF30783-C93D-4204-9703-C24C2349BA15	Esophagus Adenocarcinoma, NOS
								8373BC4A-490D-4050-BBB1-48B6472EBA34	Esophagus Adenocarcinoma, NOS
								77B3EC30-1E9F-466A-944B-6B187635EDB7	Esophagus Adenocarcinoma, NOS
								B2B8C6BC-509B-4C5B-B0FB-D07A5628FA80	Esophagus Adenocarcinoma, NOS
								1946F60D-51CA-48BA-8E94-79142594986E	Esophagus Adenocarcinoma, NOS
								CCCD0568-E267-44A2-8421-6471736BB45F	Esophagus Adenocarcinoma, NOS
								9A1C83CD-D746-4008-8ED3-22CD44951E29	Esophagus Adenocarcinoma, NOS
								608EF42F-CCD7-4AE0-8156-8775A0F09BCF	Esophagus Adenocarcinoma, NOS
								9112A681-E5A2-4DC8-B5AF-0187D3BFC1D3	Esophagus Adenocarcinoma, NOS
								48633BB8-C5F3-45AA-A406-E46C624391B4	Esophagus Adenocarcinoma, NOS
								65CEE80A-4C0A-4594-AFCA-D7088F7E5857	Esophagus Squamous Cell Carcinoma
								CD2F37A3-1912-4E64-89C5-A92A070EFD5B	Esophagus Adenocarcinoma, NOS
								AAC385C1-42F9-4D9E-B9BF-94CA13EBC5AE	Esophagus Squamous Cell Carcinoma
								3B076D77-6276-4233-867B-6A10D9E17343	Esophagus Adenocarcinoma, NOS
								89AAB5C2-9761-43BB-9053-EED42EE08969	Esophagus Squamous Cell Carcinoma
								AD2FE57C-2F44-4D57-83A9-349C8A25C414	Esophagus Adenocarcinoma, NOS
								6C7D464E-09F1-4478-837D-BFBF6ADDF062	Esophagus Adenocarcinoma, NOS
								705F9841-988A-488A-8CDD-3E00BC731E76	Esophagus Squamous Cell Carcinoma
								D86EAB07-8CB2-49A8-B00F-2986EFD7EA3B	Esophagus Squamous Cell Carcinoma
								49A3C7C6-81B0-4B14-8E13-03F6901171F8	Esophagus Adenocarcinoma, NOS
								54E45B43-363A-47A3-B43B-5E41FF589B47	Esophagus Squamous Cell Carcinoma
								D064C0BF-ECAC-4069-AEAE-4CC36BC87B38	Esophagus Squamous Cell Carcinoma
								B5FA3AC1-0691-4B6A-B6E4-765311B419ED	Esophagus Squamous Cell Carcinoma
								71F9F567-3687-4E3F-AD20-AED026A5B7E1	Esophagus Squamous Cell Carcinoma
								0B831729-8C3F-472B-9D76-50BBB4AEAB0E	Esophagus Squamous Cell Carcinoma
								C305AA0C-430C-4E3F-93AC-E611B32DB30A	Esophagus Squamous Cell Carcinoma
								FD62F7B4-A2DA-4ACD-8A85-8FA8F91521CB	Esophagus Squamous Cell Carcinoma
								2795BADD-5C70-442D-8780-8F458539D506	Esophagus Squamous Cell Carcinoma
								02D98B0A-8360-49B9-BBD6-B70C79CBD32D	Esophagus Adenocarcinoma, NOS
								DF5BD25C-D70B-4126-89CB-6C838044AE3B	Esophagus Squamous Cell Carcinoma
								91F9AD2E-E32E-43BB-AEF8-5845C398606A	Esophagus Squamous Cell Carcinoma
								484B7B36-15DC-43A9-801D-580216A8CFCB	Esophagus Squamous Cell Carcinoma
								A34C2EDA-3E24-4EB8-94CC-A5CCC819D771	Esophagus Squamous Cell Carcinoma
								CAF9B681-AF70-4A6D-A871-C45926220DE5	Esophagus Adenocarcinoma, NOS
								F4CB24E9-E7FD-4245-9C9B-7E0934038446	Esophagus Squamous Cell Carcinoma
								59DD907E-C674-46C2-BCE7-63517D5AE7A7	Esophagus Squamous Cell Carcinoma
								8D3CCAF0-872A-4EE3-B4B7-847180FEC3AB	Esophagus Squamous Cell Carcinoma
								B997D486-8C71-4926-9DF1-D2AFDAF2785F	Esophagus Squamous Cell Carcinoma
								1116DDAE-76B0-48D6-869C-BE65D8B09D35	Esophagus Squamous Cell Carcinoma
								1C2837D4-BD2C-46F7-9E22-683B82A16350	Esophagus Squamous Cell Carcinoma
								2659BB79-95F9-4476-A516-F3800536422D	Esophagus Squamous Cell Carcinoma
								FDCDF592-66F1-4027-962F-A64D49BE474F	Esophagus Adenocarcinoma, NOS
								f07070c0-fd0a-4c19-ba1e-5f06b933cd7c	Stomach  Adenocarcinoma  Diffuse Type
								6e03b415-84a1-4b91-8717-1a41edd4a255	Stomach  Adenocarcinoma  Diffuse Type
								33503349-ba7a-43f4-bcb2-ac8d16808432	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0d4d23ba-066b-4244-a3a5-922cdcaaa28f	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								12071070-d06b-4c7b-a7b1-77c84963abe6	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								151257d7-ba5b-49a0-94f5-50e28f9f46f1	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								74f859cf-6df5-4e96-ae6b-cfaab691d5fd	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								8520654b-2192-46cf-9d72-1ddb885a7c55	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								411257a4-a148-4606-8c8b-f9f1240e40ed	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								956017ba-5e4b-42c7-b93d-42316874ba93	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								7b0d69ca-13a6-46e0-839c-e433bbede3ad	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								da112431-5579-4fd1-a230-9144b359b0e9	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0d0f0200-25b8-483b-b04a-82a7833b4d30	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								117eb38b-0c62-4336-a866-fa5bd013256a	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e0916597-e489-45ef-a8f7-653a3c6591ef	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								703d3e86-32f4-44ae-bd88-c02378fc2269	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								acab1362-c8a0-448b-8533-0cd84d97959c	Stomach  Intestinal Adenocarcinoma  Tubular Type
								bec83dff-2bba-4467-92aa-5dc3db9e0eaf	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								29ed984c-74b5-4016-a649-c32e5d1d9869	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b1f261e9-c896-46d2-93c7-0290970196f8	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								66bca2f8-f269-454d-abe2-ce51a0704645	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								525a6978-2d7b-4f3a-ac68-13fe8cbbbdcc	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								d139bd0c-2eeb-44ec-856b-e77cb03d294e	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								449bf2be-b20e-4ad0-9ffb-910d83114018	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								fe9a5fcf-ff99-4dd1-8cc7-ed9f4b9d952b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								a2a71cb0-9dbd-45cc-8867-da84500515e7	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								c2f90f1a-a81b-465e-8e6b-5f02a21b27dd	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								27988b3e-cde1-4f4b-82a1-6d8ad08db8e9	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								ae633e83-303f-4a90-ac3b-34458fd9beef	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								ac2cbce4-4ba8-42c9-ac61-70f42995588e	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								75a29059-4705-4a55-bc3d-f5923b24b358	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								1ca3c5e0-32b0-4467-8ec0-ca212e35d2b3	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								05ce16d4-a999-491a-91ad-9449a57228ff	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e516880b-c4e5-43cd-a358-be0db88564d7	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								9d7305be-60a3-4ecf-856a-bd5321166590	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								af1086cc-b330-4ad6-ae1c-3b5c7259949a	Stomach  Adenocarcinoma  Diffuse Type
								daecea36-b379-46ce-8ae9-a38d22556ce2	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								43f0421c-6578-4d77-9314-ef0305c7efbd	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b14f2fcb-e869-475e-bc7c-adaf6e7e636f	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								93aebbd4-19cf-495a-b8ce-3c4a3feba152	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								df3df249-86c2-43d4-abe1-796cd7df653f	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								012e99fe-e3e8-4bb0-bb74-5b0c9992187c	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								24fed326-bdcf-4c20-a06e-7c1c3d6c9cc5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								7c831a1d-3253-4af4-a0ff-7fe9ff9af1ce	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								2f90209b-7464-4d7e-83d1-67850c263ded	Stomach  Adenocarcinoma  Diffuse Type
								44799c67-61cd-4f3e-bdbc-423e2e0fd2e8	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								1b641270-18f8-4257-821b-d26179e08e28	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								79dbe26f-930b-420a-9a2e-483e3318ca10	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								2943375e-3de8-4d42-959a-62f4192067ae	Stomach  Intestinal Adenocarcinoma  Papillary Type
								84d2d341-9e52-49ab-8c51-417208834c5d	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								69704642-ba65-4a45-8f13-95d847185045	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								ca21c55f-45a7-42af-be07-430990d5445b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								93a337ae-2bd3-4464-b38f-93dff92d3fde	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								8d3e3d09-bc61-4d5e-b8d6-5aee6cb488ae	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								bb3729fe-2fa2-47fc-afde-c14036cb84f5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								08c73b36-ede5-4136-b73b-0e964d432f89	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								9b5e5122-8ea0-432c-bc64-932bf376fffd	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0dfa566e-8fb2-4288-81ce-b4bcce674b06	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								1fb6c06e-ecf4-4b4c-bb15-df8c0a9fad8f	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b16167d0-d7d7-408e-bd0e-ce429b7f828c	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								4685dc14-1809-47ea-aba9-9bde72f6661b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								bae3c762-d894-4742-860c-43d18659b4bc	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								faed873e-a019-49c0-95a9-e9880fc23093	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								ea600f9e-cda8-4a56-b8bd-e4c01ab724ba	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								d0b78f02-2315-4c12-8b2a-ac1a110f882d	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								1d3a03a5-b757-4617-b65a-50f429a0e026	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								a89e46f4-4808-41e5-a047-acebf8744eb7	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								9f4e471e-2397-4604-b0b9-c7e0c45a5d69	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								7d915d11-729c-4577-a948-9824a8b033ce	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								889cb779-36a5-4ac6-b2cf-ea7be0cc387d	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								264b781d-122b-433d-80d4-59fbec262f22	Stomach  Adenocarcinoma  Diffuse Type
								62d32d19-ec93-4359-ab17-dc32e599b29a	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								2688d563-b052-440b-a810-873b1cd3666d	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								304aab66-4a2a-4b9f-af87-4a591d518ca1	Stomach  Adenocarcinoma  Diffuse Type
								9db0da77-797d-45d2-abdc-a586dea94ed8	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								86cbca52-29bd-4b88-a203-e8b9a072bd4b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								bdfd84fa-a8b9-46cd-b618-1267590da864	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								65ef7ccc-93e9-487a-9d4c-70ee7b9ddb17	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								632f51b3-8ef4-4b6e-a992-c81e554507ed	Stomach  Adenocarcinoma  Diffuse Type
								55a2c884-49a7-4d1a-9e0e-42f6a8b720fa	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								66c04b8f-42fd-44c5-945b-ee149bcf5dad	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b12d9857-9ae1-445b-a963-b630b27b254e	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								00781a96-4068-427c-a9c5-584d167c3dea	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								f89ac782-e0b2-4e62-844d-2e00d10a1c65	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								cae32f61-6b84-4122-af33-d8cbf410d4c8	Stomach  Adenocarcinoma  Diffuse Type
								4020b1b1-576d-4869-9ff5-552e3afb3ab5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								5b52f882-df12-42ec-a000-ff3b428ba05b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								60879669-96c2-4710-82e0-0ae9d0e5b242	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								9ebc8b58-363b-4638-8365-cfbde31a71ae	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								aef8d98f-7a1c-47f3-b5a7-a3cb299c70de	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e617b5c5-55ca-47a7-8497-854ab2844818	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								698a89d8-b4bc-48a8-80af-8414dbacf327	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e17cf9e9-0902-44e0-962b-af5492876f7d	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								fa00a6a7-d5b0-40b9-8360-d56f7db64d34	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								1a3e7ec6-c9be-4174-8799-72267eb8a7bf	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b8a1c8c7-3945-4d84-8f58-3a1a1eb10730	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								bd2a14e7-bb30-48d0-a585-7fbab767d4cf	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								ad73529b-92d7-412f-b001-920a61cec313	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								a1867d3e-a006-4e04-8118-478792bdca4b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								d59824bd-816e-4ffc-8aae-1c3eff1f800e	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0f3e2f16-1c06-4128-af45-ee73de19ea69	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								a4cd4ff1-0c0e-488d-91d9-0ac8ef3e7bfd	Stomach  Adenocarcinoma  Diffuse Type
								b604b620-cd32-4c87-ace3-00a68ab16973	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								e7a3574b-1fcf-4366-adf8-f4330d9be6f8	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								1770adbd-b424-4a88-95f0-768e5f168cee	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								03e3082e-9f8b-4689-9013-0a9954146dfa	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								ab76544e-df6a-4b05-9317-054812474d4c	Stomach  Intestinal Adenocarcinoma  Tubular Type
								d5efcde2-21e6-4071-9aa4-3c5688a45b3f	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								84d32a4c-023b-4563-9eb5-e2f47818142a	Stomach  Intestinal Adenocarcinoma  Tubular Type
								7c7cef85-024c-48c0-8644-92c034432e8f	Stomach  Intestinal Adenocarcinoma  Tubular Type
								621411ee-c1f2-4da4-b1c0-cbe0936060f6	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								337819d2-5281-4b41-9583-a476924bc837	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								98BF4652-0088-4ECC-9DFA-376438CB063C	Stomach  Adenocarcinoma  Diffuse Type
								7840D10B-4927-4731-8DDA-5F288DA21BEF	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								722A3067-2B68-442F-A333-540AB085AF32	Stomach  Intestinal Adenocarcinoma  Tubular Type
								1A1D4702-254B-4FFF-9609-182442EDB83E	Stomach Adenocarcinoma  Signet Ring Type
								10A83FCB-AD7E-48A4-A35B-6A2CBD250707	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								794FBAD5-93CE-4F9A-AEE4-37A89A6719EC	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								2A790F24-3563-4A7E-B445-56B71C70E685	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								AF59B21E-2805-4C0B-8AAD-C99B26C79EE1	Stomach  Adenocarcinoma  Diffuse Type
								30425882-20D7-47ED-9630-60977BA62C72	Stomach  Adenocarcinoma  Diffuse Type
								98C56F92-547B-473D-9FE0-2135DD5A96A8	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								CADC555E-E350-40CF-8674-3536942F3E02	Stomach  Adenocarcinoma  Diffuse Type
								E9A98A44-83F2-490C-B053-1E953EBD4E7E	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								12F63577-7FAE-4FC8-9106-39729669FDFE	Stomach  Adenocarcinoma  Diffuse Type
								8A24F1A3-4BAB-4591-9676-7BEE01D096C1	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								5A12857B-4132-4089-A7B7-E706366102DB	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								0F344863-11CC-4FAE-8386-8247DFF59DE4	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								394F390A-563C-447C-8CAE-CD7D0D7E689F	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								1EB2038D-3B4F-4DA2-964B-E629EA068215	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								5034CC97-9EEC-42FE-ACF9-CDC5A10E6547	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								DA2E8310-F54A-477E-A1FB-3F5542452780	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								FCE4E670-84BE-4592-AE6C-2B7A59218110	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0F6B5086-4BC0-4CE3-84E7-EF8EA0356D1C	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								8BB9CFD8-7012-41B1-8AA1-B6592B0C0D43	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								AA2C5D29-2764-4804-9181-5870C0956904	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0E99C305-C69F-4FB4-84D6-29CAC1CC9D39	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e42f45c6-fd00-44cf-a210-c8803da326a1	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								c991697f-77f9-412a-9bd1-0e47f7387cb2	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								eeb114ea-ecee-49ee-a162-b3b7bc0a60cd	Stomach  Intestinal Adenocarcinoma  Papillary Type
								df822ea3-e0e6-4ad3-a228-e20e3a737dac	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0aecac64-5982-4d76-8f31-958f6a00951d	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								33901c6f-9180-4a19-866c-f4ac6bec76da	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								9fa752b7-9873-4bad-9eba-8cc113705fed	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								152b9d8a-19fd-4665-b9e4-4fedfb1608cd	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								315b8de8-9842-4738-b23b-1d61722f9583	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								39acb803-f282-422e-ac56-a3a50693d203	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								6a9b3d65-ae15-4e3c-84ed-468e2196372b	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								1b98b2da-c5bf-422f-bb73-1190841fbec4	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								1a1f2197-e303-46e4-8871-53f2eb2e599a	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								5b7e557a-6aee-4077-b3cd-c96179e643ae	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								0a11e958-7ec0-493c-bd7a-eb17536facd7	Stomach  Adenocarcinoma  Diffuse Type
								efe8b105-b170-4679-8b97-c1142fefe069	Stomach  Adenocarcinoma  Diffuse Type
								8e0f35e8-cd02-4d8b-9145-597cc65eac13	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								65232ec1-4614-49bd-a3cf-d901d4f9556a	Stomach  Adenocarcinoma  Diffuse Type
								569dd23a-d0d7-409c-bc83-7340512f4195	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								b2f6d8f3-3182-4e5e-be20-8470cfb832c5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								C63AE423-C90A-4CF4-AA42-653BB8FA3E69	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								B7D2401E-7D64-4C89-A74C-16D634CAFE7F	Stomach  Adenocarcinoma  Diffuse Type
								7F056FC5-3848-4884-9644-95B8B6ED9642	Stomach  Adenocarcinoma  Diffuse Type
								A6635058-4498-463C-9BBC-1009DDA3308C	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								EC02B33C-8094-45E8-B679-B14DBE75D66F	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								028632ED-D75A-49BB-B579-41962BA502FE	Stomach  Adenocarcinoma  Diffuse Type
								23451AB5-A32B-4EEC-9B8B-FD282282A0D4	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								CDA6D95C-BE6A-485E-BA9A-C57F3EB3F99A	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								287059FB-4E69-4D6D-99D5-91EC92D35BDE	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								8ed901f6-d1bc-48b1-8514-6e8c066dc160	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								7e13fe2a-3d6e-487f-900d-f5891d986aa2	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								8e005891-d6fe-4fee-a699-971fee4e70f1	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e0e35051-fbb4-4389-aedb-87bde84d21ad	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								73d177f7-8c0a-42f6-a73d-8073c0e6142b	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								7938ca55-766e-4277-9e7e-603ffa02e570	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								6244e506-d163-4d34-8fec-6158191679cb	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								0246823f-99f4-4a46-b2b1-27f81179584c	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								8ca597db-ac21-43d1-9329-586b2eaf6982	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								37e22429-efd7-4d6a-b79a-9d8c7be2170c	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								8fccff71-3710-4cb9-b8af-61707aebe5d8	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								8cb0144b-be6b-40a1-86a2-708f96d9b615	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								5c087d33-0bd3-407d-8350-fabba6aec597	Stomach  Adenocarcinoma  Diffuse Type
								10b5774d-7079-4b4b-b824-f4745f5598c7	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								8a173d98-20a1-4c84-86c1-97818e1c665a	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								e9ea3687-0db3-47f3-8099-b4fc054397f1	Stomach  Adenocarcinoma  Diffuse Type
								c9ab5bf4-092a-420f-9b90-87deea58104c	Stomach  Adenocarcinoma  Diffuse Type
								2f1b96b2-7031-4f7d-9490-14728d320f3e	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								2ac986f6-cd6f-4ca5-90cc-7a105cfe85ea	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								306f0222-8734-44ca-9a31-86e4e62b0a1b	[Not Available]
								39dd61a8-cd17-4df3-ae17-61e239eb00bd	Stomach  Adenocarcinoma  Diffuse Type
								cad2adc1-0f41-4651-9c80-ae9961e2dc8d	[Not Available]
								442d7415-d095-4962-8296-a14390bad40a	Stomach  Adenocarcinoma  Diffuse Type
								18d0aae7-9694-4252-b8fe-16e62dadab07	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								bdc3899e-15ae-421f-829b-0e07fbd27c60	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								d5c58161-05ed-4472-81ed-176216f49929	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								263f67e0-0a28-42e6-b3f3-1ddd0d397220	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								7c9ea4fa-4cbc-4941-945a-e531e1d48304	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								588abaea-ab16-42f4-9457-5901ee791b5f	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								f35e8fef-f61d-4614-b73f-f7e1d0f21f52	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								0f00ad32-e4be-49d4-a02a-d233bfc913c5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0b4e8d73-459d-4930-9d1a-abb8d53e1844	Stomach  Adenocarcinoma  Diffuse Type
								863ad936-e2a0-471c-8d95-3a548bb94afa	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								b75eb694-8e9d-4dbb-b536-7d278c728cbb	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								1286ac33-a466-4f58-900e-f6b13fe9c25a	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								09587f1c-5c99-4102-bc49-84d50fa8d0ce	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								e62d5f3c-0a99-4932-a589-bff4fa02b1d3	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								8a9afa5e-63ce-41eb-aafb-b0b03653da5a	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								2a227292-c73d-45fc-ab68-cdfaf7aef9e7	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								2dc154aa-c478-404b-a92e-86abcc8dbdb7	Stomach  Intestinal Adenocarcinoma  Tubular Type
								a3c1f9a2-9174-48cf-81dc-2afcfdc4ba9c	Stomach  Intestinal Adenocarcinoma  Tubular Type
								8b5746f9-dbee-40bd-9141-1081960cf286	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								40e507fb-e285-424b-a98b-362e47a53e4a	Stomach  Adenocarcinoma  Diffuse Type
								e37cae9e-0656-4b89-b326-739a5cc859e0	Stomach  Adenocarcinoma  Diffuse Type
								3c573af6-6298-49de-8de0-6a7f9ce9b000	Stomach  Intestinal Adenocarcinoma  Tubular Type
								6de6b8e2-03eb-4f32-b5a3-bea73fc632b4	Stomach  Adenocarcinoma  Diffuse Type
								8bb3facb-1933-4350-9f78-74537559e000	Stomach  Adenocarcinoma  Diffuse Type
								7ef3f065-1928-4a58-b5a6-76dddb8488cf	Stomach  Adenocarcinoma  Diffuse Type
								82f9a197-040d-4c5e-89d9-d87a7aa0168d	Stomach  Adenocarcinoma  Diffuse Type
								a28e57e2-9071-4558-856a-162f7a837380	Stomach  Intestinal Adenocarcinoma  Tubular Type
								cb622cd6-776d-4681-8a67-4addd61c66b8	Stomach  Intestinal Adenocarcinoma  Papillary Type
								35496e4e-598f-4ffc-b11e-8edaeb676079	Stomach  Intestinal Adenocarcinoma  Tubular Type
								53086e20-a727-40d6-b42d-8030109b5130	Stomach  Intestinal Adenocarcinoma  Tubular Type
								6ffefb3f-3919-4fef-ab15-5d32d8f51d57	Stomach  Intestinal Adenocarcinoma  Tubular Type
								aa2934f9-8783-428e-88a6-1d04b93fa1fd	Stomach  Adenocarcinoma  Diffuse Type
								1142251c-28eb-48bf-9c54-7126d4557605	Stomach  Intestinal Adenocarcinoma  Tubular Type
								e2aaabd9-f4ba-4763-a28a-c14d46a4bc5d	Stomach  Intestinal Adenocarcinoma  Tubular Type
								1e1cdf72-820f-48a4-8fb6-89575d07cde9	Stomach  Intestinal Adenocarcinoma  Tubular Type
								8097fddc-c7b9-4192-89d7-6bcceb7465b5	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								e4344668-5a50-4dde-8eec-f4d7f01f99fd	Stomach  Intestinal Adenocarcinoma  Tubular Type
								a873475c-b2c7-4c84-bd89-5b9ee6060ce1	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								00ee3098-1b32-4e7a-81ea-993773587c41	Stomach  Intestinal Adenocarcinoma  Tubular Type
								ec8268c7-943d-451a-9036-f0196078b6c5	Stomach  Intestinal Adenocarcinoma  Tubular Type
								50279b37-828d-4e18-b41f-f793c3923e75	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								81a7f6a2-e44c-43d6-b687-01158eaf44ee	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								6D5DF29C-F2CC-43CE-8D2F-DF1DA313AD52	Stomach  Intestinal Adenocarcinoma  Tubular Type
								32B44138-64E2-4BBF-AA2F-93A3B07EC670	Stomach  Intestinal Adenocarcinoma  Tubular Type
								47FD3741-28EA-427B-9C2A-8CB39153C8AB	Stomach  Intestinal Adenocarcinoma  Tubular Type
								2B4F014F-74E9-4DDA-A1AD-ABCBD124FC44	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								17F63789-2EE3-4CD9-9273-5FA73662DB39	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								DCE21359-004F-4E8B-89E4-9E1A30FAFFDF	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								2d291fa3-f0fd-4bc1-91a6-783863714190	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								4e921bc8-df5d-4f54-93d1-d5a433f91bdb	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								FE1AD544-3189-4495-86BE-0C68F3234238	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								c3cba26e-38c2-4569-afd6-6b045e604808	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								106f3344-f699-4b2c-8f62-4ac6c948dde4	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								283776cd-32e8-4a79-bed1-1411c9d3a9e0	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								CCA74404-CD91-4A69-AE09-BFB5D9142D3D	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								1d43af3c-4279-4aed-9ced-ce3896a48ae0	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								45aacbb4-0fd8-431c-a57b-279ea57a19fc	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								d384dd0c-1947-4ffc-bee9-de38311d8249	Stomach  Adenocarcinoma  Diffuse Type
								b2e2d35c-7a79-4165-97a5-82d5c37d94f2	Stomach  Adenocarcinoma  Diffuse Type
								5956529e-f359-47c9-8222-0bd7f68c1e65	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								233cb00f-b470-47fb-ad93-2a1597655ed6	Stomach  Adenocarcinoma  Diffuse Type
								92654ad8-156c-4dd1-862e-dc339240d721	Stomach  Adenocarcinoma  Diffuse Type
								0b1663fc-6a7e-4698-951d-a8c90fdf2411	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								1a4313e7-a770-4a0e-93f1-c10562f521e3	Stomach  Intestinal Adenocarcinoma  Not Otherwise Specified (NOS)
								48ABE0EF-DFD5-4E36-A319-C4EC3165CBA1	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								0AE2B51D-AADE-4D24-8D0B-6990123D9074	Stomach Adenocarcinoma  Signet Ring Type
								6c3687fd-0761-477d-95ab-73e11d66dbd6	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								5d33895e-a959-4d5e-8850-c49c7fa49a22	Stomach  Intestinal Adenocarcinoma  Mucinous Type
								647dc54f-a51b-40ea-9d97-c93598c2af71	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								83d681f4-6c0e-45cc-a31d-0b1b31c8e073	Stomach Adenocarcinoma  Signet Ring Type
								77a8eab6-f6a1-4739-9031-75ead40d68cb	Stomach  Intestinal Adenocarcinoma  Tubular Type
								4e9d13da-438d-4b5e-8e86-67432b5e9471	Stomach  Intestinal Adenocarcinoma  Papillary Type
								e302f643-eed3-4f07-a7fb-c32ee2031a1c	Stomach  Intestinal Adenocarcinoma  Tubular Type
								edce7aac-1d8b-41cb-88e2-f50b8eeca4d8	Stomach  Adenocarcinoma  Diffuse Type
								43945bf6-6e64-4f13-9822-baf30dcd726c	Stomach  Intestinal Adenocarcinoma  Tubular Type
								c6316fb5-99a7-4adf-832f-dce44ddcd3aa	Stomach  Adenocarcinoma  Diffuse Type
								131fb0f0-e0af-431c-8462-f1925a5f4c72	Stomach  Intestinal Adenocarcinoma  Tubular Type
								3B17C66A-34AE-4B10-B907-4464D6D58383	Stomach  Adenocarcinoma  Diffuse Type
								6EEEE564-1057-42E4-BFCE-B6FE17D58D10	Stomach  Intestinal Adenocarcinoma  Tubular Type
								30A82E9F-2A2D-4A66-BDB1-26B881C21D01	Stomach  Intestinal Adenocarcinoma  Papillary Type
								16286EBD-6811-4E53-B4D5-82E1E825BC7A	Stomach  Adenocarcinoma  Diffuse Type
								8114D553-4025-4C8F-AE02-5ED01720723A	Stomach  Intestinal Adenocarcinoma  Tubular Type
								13E5A0AD-6D0C-4FB4-BF62-325F5C8ED329	Stomach  Adenocarcinoma  Diffuse Type
								AA88CE53-8439-40BD-8D93-E1D903D1C7FB	Stomach  Intestinal Adenocarcinoma  Tubular Type
								DF89B0BE-FBE1-4CA7-BEA7-AE057C3F6461	Stomach  Intestinal Adenocarcinoma  Tubular Type
								D9424178-9CB4-4BFB-9814-BA1AFE7AD954	Stomach  Intestinal Adenocarcinoma  Tubular Type
								F8A32CB0-8C3D-4BC2-8E53-561B89000BE7	Stomach  Intestinal Adenocarcinoma  Tubular Type
								C46FF2BB-BDD0-48CD-B2AE-7EF421B4B91A	Stomach  Adenocarcinoma  Diffuse Type
								7739F659-0E25-4404-AD25-EE52FBBCFD44	Stomach  Intestinal Adenocarcinoma  Tubular Type
								7997F3D8-BFF6-4808-8DF0-D88D344D7C9E	Stomach  Intestinal Adenocarcinoma  Tubular Type
								B93221BA-2076-4DD4-9931-179571B7E16C	Stomach  Adenocarcinoma  Diffuse Type
								75915370-0895-4933-8875-8A31C65C3D50	Stomach  Adenocarcinoma  Diffuse Type
								41327AAC-1228-435A-8DE5-B26886BDBCE0	Stomach  Adenocarcinoma  Diffuse Type
								FE45EAEE-FF0C-4D52-B18C-04DC2D64D3B4	Stomach  Adenocarcinoma  Diffuse Type
								7B594344-573D-46A0-B735-6727F6F15585	Stomach  Intestinal Adenocarcinoma  Tubular Type
								39FA91BD-DFFD-4B2C-AF7D-7FD914730B08	Stomach  Adenocarcinoma  Diffuse Type
								00FDF96D-C48A-47E0-B91D-576C52E6819C	Stomach  Intestinal Adenocarcinoma  Tubular Type
								D262FDF0-30BF-4618-9342-DE5A2A053828	Stomach  Intestinal Adenocarcinoma  Papillary Type
								416FD59B-A4F8-4610-9F8B-473CB5C04A4A	Stomach  Intestinal Adenocarcinoma  Tubular Type
								98760B3E-3075-4444-9A01-8F11935069D4	Stomach  Intestinal Adenocarcinoma  Tubular Type
								D807A1C0-F586-4F13-B8EF-3F79AB3D4B70	Stomach  Intestinal Adenocarcinoma  Tubular Type
								87c217d4-66f4-46ca-8244-7856ce658fd3	Stomach  Adenocarcinoma  Diffuse Type
								8f9be4e3-a1e5-4fa3-8d31-a53de781ed97	Stomach  Adenocarcinoma  Diffuse Type
								ba3911e6-2208-4a04-8891-b5cba9fda6bf	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								d5550422-dc4b-4ebf-acb5-aa5243fff49f	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)
								faf5f1e3-9f7d-4a9d-a005-6cad5e4becf7	Stomach  Adenocarcinoma  Not Otherwise Specified (NOS)"""
		reader = [i.strip() for i in harcoded_histology.splitlines()]
		reader = [i.split('\t') for i in reader]
		reader = [(i[0].lower(), i[1]) for i in reader]
		reader = dict(reader)
	return reader

def _file_api_basic_info(info):
	basic_info = dict()
	if 'data' in response :
		case = response['data']['cases'][0]
	
		#Tissue Type
		tissue_type = None
		if 'samples' in case:
			tissue_id = case['samples'][0].get('sample_type_id', None)#Defaults to 'control' for tissue type
			if isinstance(tissue_id, str):
				tissue_type = ('tumor' if int(tissue_id) < 10 else ('normal' if int(tissue_id) < 20 else 'control'))
		
		#Read Groups
		try:
			read_groups = response['data']['analysis']['metadata']['read_groups'][0]
		except:
			read_groups = None

		try:
			index_file = response['data']['index_files'][0]
		except:
			index_file = None
	else:
		tissue_id = None
		tissue_type = None
		read_groups = None
		index_file = None

	info['basic_info'] = {
		#---------------- Basic file info ----------------
		'file_size': response['data']['file_size'],
		'file_id': response['data']['file_id'],
		'file_name': response['data']['file_name'],
		#------------------- sample info -----------------
		'barcode': (case['samples'][0]['submitter_id'] if 'samples' in case.keys() else ""),
		'tissue_code': tissue_id,
		'tissue_type': tissue_type,
		'genome_build': "GRCh38.p0",
		'genome_name': "GRCh38.d1.vd1",
		#---------------- Basic case info ----------------
		'case_id': case['case_id'],
		'patientID': case['submitter_id'],
		#--------------- Probe Targets info --------------
		'exome targets link': read_groups.get('target_capture_kit_target_region', ""),
		'capture kit catalog number': read_groups.get('target_capture_kit_catalog_number', ""),
		'capture kit name': read_groups.get('target_capture_kit_name', "").replace('\n', ''),
		#Index Files
		'index_id': index_file.get('file_id', ""),
        'index_name': index_file.get('file_name', "")
	}
	return info

def file_api(file_id, parameters = None, use_local = False, legacy = False):
	""" Requests information about a file from the GDC API
		Parameters
		----------
			file_id: string
				The uuid for the file
			legacy: bool; default False
				Whether to use to GDC legacy archive
			fields: list<string>
				A list of fields to request. Only pass top-level variables,
				all subvariables will be returned.
				Available top-level fields:
					associated_entities 4
					downstream_analyses 28
					metadata_files      16
					analysis            65
					annotations         17
					index_files         17
					archive             15
					center               6
					cases              300 (Currently not supported)
		Returns
		-------
			* 'access': string
			* 'acl':    list<string>
			* 'analysis': dict<>
			* 'associated_entities': list<dict>
			* 'basic_info': dict<>
				* 'aligned_reads': list<dict>
					* 'barcode': string [Sample Barcode]
					* 'category': string [Histology]
					* 'experimental_strategy': string
					* 'filename'
					* 'id'
					* 'md5'
					* 'patient': string [Patient Barcode]
					* 'size': int
					* 'state': string
					* 'tissue': string {'Primary Tumor', 'Solid Tissue Normal', 'Blood Derived Normal', 'Metastatic'}
				* 'barcode': string [the sample barcode]
				* 'capture kit catalog number': string
				* 'capture kit name': string
				* 'case_id': string [Case UUID]
				* 'exome targets link': string [""]
				* 'file_id': string [file UUID]
				* 'file_name': string
				* 'file_size': int
				* 'genome_build': string
				* 'reference_name': string
				* 'index_id': string
				* 'index_md5sum': string
				* 'patient_id': string [Patient Barcode]
				* 'sample_type': string {'Primary Tumor', 'Solid Tissue Normal', 'Blood Derived Normal', 'Metastatic'}
			* 'cases': list<dict>
			* 'created_datetime': string
			* 'data_category': string {'Raw Sequencing Data'}
			* 'data_format': string {'BAM'}
			* 'data_type': string   {'Aligned Reads'}
			* 'downstream_analyses': list<dict<>>
			* 'experimental_strategy': string ['WXS']
			* 'file_id': string
			* 'file_name': string
			* 'file_size': int
			* 'file_state': string
			* 'index_files': list<dict<>>
			* 'md5sum': string [UUID]
			* 'platform': string
			* 'state': string
			* 'submitter_id': string [UUID]
			* 'type': string ['aligned_reads']
			* 'updated_datetime': string
	"""
	file_id = file_id.lower()
	if legacy:
		endpoint = "legacy/files"
	else: endpoint = 'files'
	if not parameters:
		all_expand_fields = ['analysis', 'analysis.metadata.read_groups', 'analysis.input_files', 
							 'annotations', 'archive', 'associated_entities', 'cases', 'cases.annotations',
							 'cases.demographic', 'cases.diagnoses', 'cases.family_histories', 'cases.project',
							 'cases.samples', 'cases.samples', 'downstream_analyses.output_files', 'center']
		all_expand_fields = ','.join(all_expand_fields)
		parameters = {
			'expand': all_expand_fields
		}

	#----------------------- Send the request to the GDC Server --------------------------
	if not use_local:
		try:
			response = _request(endpoint, file_id, parameters)
		except:
			use_local = True
	if use_local:
		with open(local_file_api_filename, 'r') as file1:
			local_api = json.loads(file1.read())
		response = local_api.get(file_id)
		
	#--------Include a key to easily access some important information on the file--------- 
	#pprint(response)
	case = response['cases'][0]
	if 'samples' in case.keys():
		sample = case['samples'][0]
		sample_type = sample['sample_type']
	
	try:
		read_groups = response['analysis']['metadata']['read_groups'][0]
	except:
		read_groups = {}

	try:
		index_file = response['index_files'][0]
	except:
		index_file = {}

	file_info = {
		#---------------- Basic file info ----------------
		'file_size': response['file_size'],
		'file_id': response['file_id'],
		'file_name': response['file_name'],
		#------------------- sample info -----------------
		'barcode': sample['submitter_id'],
		'sample_type': sample_type,
		'genome_build': "GRCh38.p0",
		'genome_name': "GRCh38.d1.vd1",
		#---------------- Basic case info ----------------
		'case_id': case['case_id'],
		'patient_id': case['submitter_id'],
		#--------------- Probe Targets info --------------
		'exome targets link': read_groups.get('target_capture_kit_target_region', None),
		'exome file': _get_exome_targets(read_groups.get('target_capture_kit_catalog_number', "")),
		'capture kit catalog number': read_groups.get('target_capture_kit_catalog_number', ""),
		'capture kit name': read_groups.get('target_capture_kit_name', "").replace('\n', ''),
		#Index Files
		'index_id': index_file.get('file_id'),
        'index_name': index_file.get('file_name'),
        'index_md5sum': index_file.get('md5sum')
	}
	#response = response['data']
	response['basic_info'] = file_info
	return response

def case_api(case_id, parameters = None, use_local = False):
	""" Retrieves information from the GDC API concerning a specific case
		Parameters
		----------
			case_id: string
				The uuid linked to the case
			fields: list<string>; default None
				Any fields to retrieve. Accepts upper-level arguments and
				finds all lower-level arguments automatically.
			expand: 'Full', list<string>; default 'Full'
				A list of sections to expand on. 'Full' retrieves expanded sections
				for all headers that are compatible with the service
		Returns
		----------
			case_info: dict<>
				The information for the case.
				Top-level keys:
				* case_info: Basic info concerning the case
				* 'aliquot_ids':           list<string>
				* 'analyte_ids':           list<string>
				* 'basic_info':            dict<>
					* 'aligned_reads': list<>
					* 'case_id': string
					* 'file_count': int
					* 'file_size':  int
					* 'histology':  string
				* 'case_id':               string
				* 'created_datetime':      string, None
				* 'demographic':           dict<>
					* 'created_datetime': None
					* 'demographic_id':   string
					* 'ethnicity':        string
					* 'gender':           string ['female', 'male']
					* 'race':             string
					* 'state':            None
					* 'submitter_id':     string
					* 'updated_datetime': string [ISO datetime]
					* 'year_of_birth':    int
					* 'year_of_death':    int, None
				* 'diagnoses':             list<dict<>>
				* 'exposures':             list<dict<>>
				* 'files':                 list<dict<>>
					* 'access': string ['open', 'controlled']
					* 'acl': list<string>
					* 'created_datetime':  string
					* 'data_category':     string
					* 'data_format':       string
					* 'data_type':         string
					* 'experimental_strategy': string
					* 'file_id': string
					* 'file_name': string
					* 'file_size': int
					* 'file_state': string
					* 'md5sum': string
					* 'platform': string
					* 'state': string
					* 'submitter_id': string
					* 'type': string
					* 'updated_datetime': string
				* 'portion_ids':           list<string>
				* 'sample_ids':            list<string>
				* 'samples':               list<dict<>>
					* 'composition': 		None
					* 'created_datetime': 	None
					* 'current_weight': 	None
					* 'days_to_collection': None
					* 'freezing_method': 	None
					* 'initial_weight':		None
					* 'intermediate_dimension': string
					* 'is_ffpe': 			bool
					* 'longest_dimension':  string
					* 'oct_embedded': 		None
					* 'pathology_report_uuid': string
					* 'preservation_method':None
					* 'sample_id':			string
					* 'sample_type':		string ['Primary Tumor']
					* 'sample_type_id':		string [INT]
					* 'shortest_dimension': string
					* 'submitter_id': 		string [TCGA SAMPLE BARCODE]
				* 'slide_ids':             list<string>
				* 'state': None
				* 'submitter_aliquote_ids':list<string>
				* 'submitter_analyte_ids': list<string>
				* 'submitter_id': string
				* 'submitter_portion_ids': list<string>
				* 'submitter_sample_ids':  list<string>
				* 'submitter_slide_ids':   list<string>
				* 'summary': dict<>
					* 'file_count': int
					* 'file_size':  int [bytes]
				* 'updated_datetime': string [ISO Datetime]
	"""
	case_id = case_id.lower()
	if not parameters:
		parameters = {
			'expand': ','.join(['samples', 'files', 'files.cases.samples','summary', 'exposures', 'diagnoses', 'demographic'])
		}
	#----------------------- Send the request to the GDC Server --------------------------
	if not use_local:
		try:
			response = _request('cases', case_id, parameters)
		except:
			use_local = True
	if use_local:
		with open(local_case_api_filename, 'r') as file1:
			local_api = json.loads(file1.read())
			response = local_api.get(case_id)
	cancer_type = histology.get(case_id)

	case_files = generate_manifest_file([response])

	case_info = {
		'case_id': response['case_id'],
		'file_count': response['summary']['file_count'],
		'file_size': response['summary']['file_size'],
		'patient_barcode': response['submitter_id'],
		'aligned_reads': case_files,
		'histology': cancer_type
	}
	response['basic_info'] = case_info
	return response

def _sample_list_from_ids(case_ids, options):
	bam_folder = options['bam folder']
	exome_targets_folder = options['targets folder']

	sample_list = list()
	exome_targets_table = list()
	for index, case_id in enumerate(case_ids):
		print("{0} of {1}: {2}".format(index+1, len(case_ids), case_id))
		case_data = case_api(case_id)
		genome_files = case_data['basic_info']['aligned_reads']

		N = [f for f in genome_files if f['tissue'] == 'Solid Tissue Normal']
		if len(N) == 0:
			N = [f for f in genome_files if f['tissue'] == 'Blood Derived Normal']

		T = [f for f in genome_files if f['tissue'] == 'Primary Tumor']

		if len(N) != 1 or len(T) != 1:
			print(case_id, " does not have the required genome files:")
			print("Available: ", ", ".join([i['tissue'] for i in genome_files]))
			print("\tNormal: ", type(N), len(N))
			print("\tTumor: ", type(T), len(T))
			continue
		N = N[0]
		T = T[0]
		#exome_targets = N['basic_info']['exome targets link']
		#exome_targets = _map_target_probe_files(exome_targets)
		exome_targets = file_api(T['id'])['basic_info']
		#exome_targets = _get_exome_targets(exome_targets)
		#print("Exome Targets:")
		#print("\tBarcode: ", case_data['basic_info']['patient_barcode'])
		#print("\tCapture Kit Name: ", exome_targets['capture kit name'])
		#print("\tCapture Kit Catalog Number: ", exome_targets['capture kit catalog number'])
		#print("\tExome Targets Link: ", exome_targets['exome targets link'])
		_et = {
			'Patient Barcode': case_data['basic_info']['patient_barcode'],
			'Capture Kit Name': exome_targets['capture kit name'],
			'Capture Kit Catalog Number': ','.join([i for i in exome_targets['capture kit catalog number'].split('|') if i != 'NA']),
			'Exome Targets Link': exome_targets['exome targets link'],

		}
		exome_targets_table.append(_et)
		#exome_targets = exome_targets['analysis']['metadata']['read_groups'][0]['library_name']
		#print("{0} of {1}".format(index+1, len(case_ids)), case_data['submitter_id'], exome_targets)
		
		exome_targets = "/home/upmc/Documents/Reference/Capture_Targets/SeqCap_EZ_Exome_v3_GRCh38_UCSCBrowser_capture_targets.bed"
		sample = {
			'CaseID':  case_data['case_id'],
			'SampleID':  T['barcode'],
			'SampleUUID':T['id'],
			'NormalID':  N['barcode'],
			'NormalUUID':N['id'],
			'PatientID': case_data['submitter_id'],
			'TumorBAM':  bam_folder + '/'.join([T['id'], T['filename']]), #designed for linux
			'NormalBAM': bam_folder + '/'.join([N['id'], N['filename']]),
			'ExomeTargets': exome_targets,
			'Genome Type': T['experimental_strategy'],
			'Histology': case_data['basic_info']['histology']
		}
		sample_list.append(sample)
		sample_list = sorted(sample_list, key = lambda s: s['PatientID'])

	with open('C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Data\\exome_targets_table.tsv', 'w', newline = "") as file1:
		writer = csv.DictWriter(file1, delimiter = "\t", fieldnames = sorted(exome_targets_table[0].keys()))
		writer.writeheader()
		writer.writerows(exome_targets_table)
	return sample_list

def _sample_list_from_folder(options):
	bam_folder = options['bam folder']
	exome_targets = options['targets folder']
	sample_list = list()
	#------------ Get data for each file and organize by patient and tissue source site ------------
	patients = dict()
	for file_id in os.listdir(bam_folder): #Basically a list of ids for the genome files in that folder
		if '-' not in file_id: continue #skip any non-genome folders

		file_data = file_api(file_id)
		patient_barcode = file_data['basic_info']['case_id']
		tissue_source = file_data['basic_info']['sample_type']

		if patient_barcode not in patients:
			patients[patient_barcode] = dict()
			case_id = file_data['basic_info']['case_id']
			patients[patient_barcode]['case data'] = case_api(case_id)

		patients[patient_barcode][tissue_source] = file_data
	#-------------------------------- Generate the sample table -------------------------------------
	for patient_barcode, patient_reads in sorted(patients.items()):
		case_data = patient_reads['case data']
		
		T = patient_reads.get('Primary Tumor')

		if 'Blood Derived Normal' in patient_reads.keys(): N = patient_reads['Blood Derived Normal']
		elif 'Solid Tissue Normal' in patient_reads.keys(): N = patient_reads['Solid Tissue Normal']
		else: N = None

		if N is None or T is None:
			print(patient_barcode, " was missing one or more required files (Available: {0}).".format(', '.join(sorted(patient_reads.keys()))))
			continue

		row = {
			'CaseID':  case_data['case_id'],
			'SampleID':  T['basic_info']['barcode'],
			'SampleUUID':T['basic_info']['case_id'],
			'NormalID':  N['basic_info']['barcode'],
			'NormalUUID':N['basic_info']['case_id'],
			'PatientID': case_data['submitter_id'],
			'TumorBAM':  bam_folder + '/'.join([T['basic_info']['file_id'], T['basic_info']['file_name']]), #designed for linux
			'NormalBAM': bam_folder + '/'.join([N['basic_info']['file_id'], N['basic_info']['file_name']]),
			'ExomeTargets': exome_targets,
			'Genome Type': "WXS",
			'Histology': case_data['basic_info']['histology']
		}
		if row['Histology'] == 'Esophagus Adenocarcinoma, NOS':
			sample_list.append(row)
	return sample_list

def generate_sample_list(io, filename = None):
	""" Generates a list of samples to use with the genome pipeline.
		Parameters
		----------
			io: string, list<string>
				Either the path to a folder with genome files, or a list of case ids.
		Required Fields
		---------------
			* 'CaseID': 
			* 'SampleID'
			* 'NormalID'
			* 'PatientID':
			* 'NormalBAM':
			* 'TumorBAM':
			* 'Targets':
			* 'Type': {'WXS', 'WGS'}
			* 'Notes': "" [OPTIONAL]
	"""
	if isinstance(io, str):
		options = {
			'bam folder': io,
			'targets folder': ""
		}
		sample_list = _sample_list_from_folder(options)
	else:
		options = {
			'bam folder': "",
			'targets folder': ""
		}
		sample_list = _sample_list_from_ids(io, options)

	if filename is not None and len(sample_list) > 0:
		with open(filename, 'w', newline = "") as file1:
			writer = csv.DictWriter(file1, fieldnames = sorted(sample_list[0].keys()), delimiter = '\t')
			writer.writeheader()
			writer.writerows(sample_list)

	return sample_list

def generate_manifest_file(case_ids, filename = None, debug = False):
	""" Generates a manifest file for all genome files attached to the provided case_ids.
		Parameters
		----------
			case_id: list<string, dict<>>
				A list of either case_ids or case_api responses.
		Returns
		-------
			manifest: list<dict>
				A list of the relevant fields used to generate a manifest file.
				* 'id': 		string
				* 'filename': 	string
				* 'md5': 		string
				* 'size': 		int
				* 'state': 		string
				* 'patient': 	string
				* 'category': 	string
				* 'tissue': 	string
				* 'barcode': 	string
	"""
	manifest = list()
	for index, case_id in enumerate(case_ids):
		if debug:
			print("{0} of {1}".format(index+1, len(case_ids)))
		if isinstance(case_id, str):
			case_info = case_api(case_id)
		else:
			case_info = case_id

		aligned_reads = [f for f in case_info['files'] if (f['data_category'] == 'Raw Sequencing Data' and f['experimental_strategy'] in ['WXS', 'WES', 'WGS'])]

		for read in aligned_reads:
			sample = read['cases'][0]['samples'][0]
			barcode = sample['submitter_id']
			patient_code = '-'.join(barcode.split('-')[:3])
			category = histology.get(case_info['case_id'])
			fields = {
				'id': read['file_id'],
				'filename': read['file_name'],
				'md5': read['md5sum'],
				'size': read['file_size'],
				'state': read['file_state'],
				'patient': patient_code,
				'category': category,
				'tissue': sample['sample_type'],
				'barcode': barcode,
				'experimental_strategy': read['experimental_strategy']
			}
			manifest.append(fields)
	if filename:
		with open(filename, 'w', newline = "") as file1:
			fieldnames = ["id",	"filename",	"md5",	"size",	"state"]
			_excluded = set(manifest[0].keys()) - set(fieldnames)
			fieldnames += sorted(_excluded)
			writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = fieldnames)
			writer.writeheader()
			writer.writerows(manifest)
	return manifest

def generate_local_api(filename):
	""" Saves a local copy of the relevant files"""
	case_ids = [i.strip() for i in histology.split('\n')]
	case_ids = [i.split('\n')[0].lower() for i in case_ids]

	for case_id in case_ids:
		pass

	return filename

histology = _get_histology()
if __name__ == "_main__":
	case_id  = "6969fe5a-5993-48e5-95c5-c5c7d3d08205"
	#file_id  = "2d4f1ce4-4613-403a-90ec-fd6a551b6487"
	file_id  = "0496553c-c68b-42a5-8a54-29f20b8f7c44"
	index_id = "4570dd4d-9234-4d37-b8d0-66d37594e3f1"
	stomach_case_id = "00781a96-4068-427c-a9c5-584d167c3dea"

	pprint(file_api(file_id))
	#pprint(case_api(case_id))
else:
	import pandas
	data = pandas.read_excel("C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Clinical Data\\nationwidechildrens.org_clinical_patient_esca.xlsx")
	all_esca_case_ids = sorted(i.lower() for i in data['bcr_patient_uuid'].values if '-' in i)

	data = pandas.read_excel("C:\\Users\\Deitrickc\\Documents\\UPMC Files\\Projects\\Genome Instability Project\\Clinical Data\\20140110_STAD_Clinical_Data_Blacklisted_Cases_Removed.xlsx", skiprows = 1)
	all_stad_case_ids = sorted(i.lower() for i in data['bcr_patient_uuid'].values if '-' in i)


	filename = "G:\\Pipeline Files\\Combined Pipeline\\Full_ESCA_sample_list.tsv"
	generate_sample_list(all_esca_case_ids, filename = filename)

	filename = "G:\\Pipeline Files\\Combined Pipeline\\Full_STAD_sample_list.tsv"
	generate_sample_list(all_stad_case_ids, filename = filename)

print("Finished!")



