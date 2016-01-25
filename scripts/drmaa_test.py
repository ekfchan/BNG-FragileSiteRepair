import drmaa
from string import Template
import time

cSession = drmaa.Session()
cSession.initialize()

jTemplate=cSession.createJobTemplate()
jTemplate.remoteCommand="drmaa_test.sh"
jTemplate.joinFiles=True
jTemplate.outputPath=":."
jTemplate.jobName="drmaa_test"
jTemplate.nativeSpecification=""
#jTemplate.nativeSpecification="-pe openmp 8"

print("Starting job")
jobID=cSession.runJob(jTemplate)
print("Job started")

i=0
while not cSession.jobStatus(jobID) in [drmaa.JobState.DONE, drmaa.JobState.FAILED]:
	print("Waiting..%s" % i)
	time.sleep(1)
	i+=1
	
print("Job complete")



