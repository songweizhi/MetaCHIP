
def read_in_job_script_header(job_script_header_example):

    job_script_header = ''
    for each in open(job_script_header_example):
        job_script_header += each

    job_script_header += '\n'

    return job_script_header


job_script_header_example = 'input_file_examples/blastn_job_script_header_demo.sh'
job_script_header = read_in_job_script_header(job_script_header_example)
print(job_script_header)
print('haha')



