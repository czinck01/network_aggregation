from __future__ import division
import sys
import os

script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
lib_path = os.path.abspath(os.path.join(script_path, ".."))
sys.path.insert(0, lib_path)

from funcassociate.client import FuncassociateClient

c = FuncassociateClient()
response = c.functionate(query=("YBL071W-A", "YCL055W", "YCR094W",
                                "YDL127W", "YFL029C", "YGR271C-A",
                                "YHR099W", "YJR066W", "YKL203C",
                                "YNL289W"),
                         species="Saccharomyces cerevisiae",
                         namespace="sgd_systematic")

print "OVERREPRESENTED ATTRIBUTES"

headers = ("N", "X", "LOD", "P", "P_adj", "attrib ID", "attrib name")
print "\t".join(headers)

info = response["request_info"]
reps = info["reps"]
below_detection_threshhold = "< %f" % (1/reps)

for row in response["over"]:
    row.pop(1)
    if row[4] is 0:
        row[4] = below_detection_threshhold
    print "\t".join(map(str, row))

print "\nREQUEST INFO"
info = response["request_info"]
for k in info.keys():
    print "%s: %s" % (k, info[k])
