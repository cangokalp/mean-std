# NETGEN-DEG instances (n = 4096, m ranges from 2*n to n*n)

#         seed   #         n   src   trg         m    costs   supply                 caps
echo "1   1      4096    64    64      8192  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_01a.min.param
echo "1   2      4096    64    64     16384  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_02a.min.param
echo "1   3      4096    64    64     32768  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_03a.min.param
echo "1   4      4096    64    64     65536  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_04a.min.param
echo "1   5      4096    64    64    131072  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_05a.min.param
echo "1   6      4096    64    64    262144  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_06a.min.param
echo "1   7      4096    64    64    524288  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_07a.min.param
echo "1   8      4096    64    64   1048576  1 10000    64000  0 0 100 100  1 1000" > netgen_deg_08a.min.param


#    seed   #      n    src   trg         m    hosts   supply                 hfps
echo "1  10      1024    32    32      32768  1 10000    32000  0 0 100 100  1 1000" > netgen_32_10a.min.param
echo "1  11      2048    45    45      65536  1 10000    45000  0 0 100 100  1 1000" > netgen_32_11a.min.param
echo "1  12      4096    64    64     131072  1 10000    64000  0 0 100 100  1 1000" > netgen_32_12a.min.param
echo "1  13      8192    91    91     262144  1 10000    91000  0 0 100 100  1 1000" > netgen_32_13a.min.param
echo "1  14     16384   128   128     524288  1 10000   128000  0 0 100 100  1 1000" > netgen_32_14a.min.param
echo "1  15     32768   181   181    1048576  1 10000   181000  0 0 100 100  1 1000" > netgen_32_15a.min.param


#    seed   #      n    src   trg         m    hosts   supply                 hfps
echo "1  10      1024    32    32      32768  1 100    320  0 0 100 100  1 1000" > netgen_lo_32_10a.min.param
echo "1  11      2048    45    45      65536  1 100    450  0 0 100 100  1 1000" > netgen_lo_32_11a.min.param
echo "1  12      4096    64    64      131072  1 100    640  0 0 100 100  1 1000" > netgen_lo_32_12a.min.param
echo "1  13      8192    91    91     262144  1 100    910  0 0 100 100  1 1000" > netgen_lo_32_13a.min.param
echo "1  14     16384   128   128     524288  1 100   1280  0 0 100 100  1 1000" > netgen_lo_32_14a.min.param
echo "1  15     32768   181   181     1048576  1 100   1810  0 0 100 100  1 1000" > netgen_lo_32_15a.min.param
