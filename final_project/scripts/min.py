#!/usr/bin/env python
# coding:utf-8


import glob
import os

#----------------------------------------------------------------------


def analysis_txt(txt):
    """"""
    f1 = open(txt, 'r')
    for x in f1:
        f11 = x.split('\t')
        ana[f11[0] + f11[1]] = f11
    return ana.copy()

#----------------------------------------------------------------------
#----------------------------------------------------------------------


def cutoff(a):
    """"""
    f2 = open(a, 'r')
    for f3 in f2:
        f3 = f3.strip('\n')
        f33 = f3.split(',')
        cut[f33[0]] = f33
    return cut
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

def comment(score):
    if score > 5:
        return "Mutation confirmed"
    elif score > 0:
        return "Ambiguous sequencing result"
    else:
        return "Mutation untrusted"

#-----------------------------------------------------------------------
if __name__ == '__main__':

    txtfiles = glob.glob("scripts/txt/*")
    for txt in txtfiles:
        if txt.endswith('txt'):
            f = os.path.basename(txt).split('.')[0]
            ana = {}
            ana_txt = {}
            ana_txt = analysis_txt(txt)
            ana.clear()
            vcffile = 'scripts/vcf/' + f.upper() + '_afterinter.vcf'
            ana = {}
            ana_vcf = {}
            ana_vcf = analysis_txt(vcffile)
            ana.clear()
            result = []
            for vcf_x in ana_vcf:
                if vcf_x in ana_txt.keys():
                    result.append(ana_txt[vcf_x])
            result_news = []
            for r in result:
                for rr in range(3, 9):
                    if int(r[rr]) > 0:
                        result_news.append(r)

            l = {}
            cut = {}
            cutcsv = cutoff('scripts/cutoff/' + f.upper() + 'cutoff.final.csv')

            ff = open(vcffile + 'new.vcf', 'w')
            for result_new in result_news:

                l['LMNA'] = int(result_new[3])
                l['TNNT2'] = int(result_new[4])
                l['SCNSA'] = int(result_new[5])
                l['MYBPC3'] = int(result_new[6])
                l['MYH6'] = int(result_new[7])
                l['MYH7'] = int(result_new[8])
                if l['LMNA'] > 0:
                    res = 'LMNA'
                    resv = l['LMNA']
                    resvv=result_new[2]
                elif l['TNNT2'] > 0:
                    res = 'TNNT2'
                    resv = l['TNNT2']
                    resvv=result_new[2]

                elif l['SCNSA'] > 0:
                    resv = l['SCNSA']
                    res = 'SCNSA'
                    resvv=result_new[2]

                elif l['MYBPC3'] > 0:
                    resv = l['MYBPC3']
                    res = 'MYBPC3'
                    resvv=result_new[2]

                elif l['MYH6'] > 0:
                    resv = l['MYH6']
                    res = 'MYH6'
                    resvv=result_new[2]

                elif l['MYH7'] > 0:
                    resv = l['MYH7']
                    res = 'MYH7'
                    resvv=result_new[2]

                # print f,
                # print 'resv:',resv,
                # print 'cutcsv[res]:',res,len(cutcsv[res])
                if len(cutcsv[res]) > resv:
                    if cutcsv[res][resv] != 'NULL':
                        if cutcsv[res][resv] == '':
                            cutcsv[res][resv] = 0
                        # print result_new,int(cutcsv[res][resv])-resv
                        # print
                        # result_new[0]+result_new[1],int(cutcsv[res][resv])-resv
                        k = result_new[0] + result_new[1]
                        score = int(resvv) - int(cutcsv[res][resv])
                        last =  comment(score)

                        if k in ana_vcf.keys():
                            # print ana_vcf[k],last
                            # aaa=ana_vcf[k][-1]
                            ana_vcf[k][-1] = str(ana_vcf[k][-1]).strip('\n')
                            ana_vcf[k].append(last)
                            # print ana_vcf[k]
                            # s=ana_vcf[k].join('\n')

                            for x in ana_vcf[k]:
                                ff.write(str(x))
                                ff.write('\t')
                            ff.write('\n')

            ff.close()
