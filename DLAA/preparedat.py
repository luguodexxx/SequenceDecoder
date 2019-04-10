#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2019/4/10 5:32 PM
__author__ = 'Zhou Ran'

"""
This is for prepare the training data
"""

# TODO, add kmer?

import re
import numpy as np
import gzip
import pysam
import pandas as pd
import random
from collections import defaultdict

random.seed(9)


class FormatError(Exception):
    pass


class AlphaToMtx:

    def __init__(self,
                 seqinfo,
                 mode):
        """
        Convert the sequence into dataframe or assay, dataframe to CNN and assay for LSTM
        :param seqlist:
        :param mode:
        """
        if mode == 'cnn':
            self.__seqlist = seqinfo['sequence']
            self.__label = seqinfo['target']
            self.__seqlen = len(self.__seqlist[0])
            self.__trainlen = len(self.__seqlist)

            self.train_mtx = np.zeros((self.__trainlen,
                                       self.__seqlen, 4))
            self.label_mtx = np.zeros((self.__trainlen,
                                       len(set(self.__label))))

            for i in range(self.__trainlen):
                self.train_mtx[i, :] = self.__cnn_mtx(self.__seqlist[i])
                print(self.label_mtx)
                self.label_mtx[i, int(self.__label[i])] = 1

            self.label_mtx.astype('float32')
            self.train_mtx.astype('float32')

        elif mode == 'rnn':
            self.train_mtx = seqinfo['sequence'].apply(
                lambda x: [int(self.__rnn_lst(e)) for e in x])
            self.label_mtx = np.array(seqinfo['target'])
        else:
            raise TypeError("only concern for cnn of rnn mode.")

    @staticmethod
    def __cnn_mtx(seq):
        """
        convert the dna alpha into one-hot encoding format
        :param seq: str
        :return: numpy.ndarray
        """
        seq = seq.upper()

        alpha = ['[aA]', '[cC]', '[gG]', '[tT]']

        for i in range(len(alpha)):
            index_ = [base.start() for base in re.finditer(alpha[i], seq)]
            tem = np.zeros(len(seq), dtype=np.int)
            tem[index_] = 1

            if i == 0:
                temass = np.zeros((len(seq), 4))
            temass[..., i] = tem

        return temass

    @staticmethod
    def __rnn_lst(letter):
        letter = letter.upper()
        alpha = 'ATGC'
        return next((i for i, _letter in enumerate(
            alpha) if _letter == letter), None)


class SeqCrawler:
    """
    1) bed region and point
    2) fasta file
    """
    global __SEQDIC

    def __init__(self,
                 site_lst,
                 left,
                 right,
                 label=None,
                 genome_fasta=None):
        """

        :param sitelst:
        :param fafaidx:
        """
        self.__sitelst = site_lst
        self.__left = int(left)
        self.__right = int(right)

        self.__checkbed(site_lst[0])
        self.__fa = pysam.FastaFile(genome_fasta)

        self.mtx = self.__getfa()
        if label != None:
            self.mtx['target'] = label
        else:
            self.mtx['target'] = 'NaN'

    def __getfa(self):
        """
        fetch the sequence from fasta
        :return:
        """
        __res = []
        for _site in self.__sitelst:
            chrom, s, e, _, _, strand = _site.strip().split('\t')
            seq = self.__fa.fetch(chrom, int(s) - self.__left, int(s) + self.__right)

            if "N" in seq: continue

            if strand == '-':
                seq = self.revseq(seq)[::-1]

            __res.append(seq)

        return pd.DataFrame(__res, columns=['sequence'])

    @classmethod
    def readbed(cls, file, left, right, label=None, genome=None, sample=None):
        """
        read bed file and return a list
        :param file:
        :return:
        """
        __res = []
        fh = open(file)
        for line in fh.readlines():
            chrom, s, e, _, _, strand = line.strip().split('\t')
            if len(chrom) > 6: continue
            __res.append(line)
        fh.close()
        if sample:
            __res = sample(__res, sample)
        return cls(__res, left, right, label, genome)

    @classmethod
    def __checkbed(cls, info):

        _tocheck = info.strip().split('\t')

        if len(_tocheck) != 6:
            raise FormatError("File don't look like bed6 format!")
        elif abs(int(_tocheck[2]) - int(_tocheck[1])) != 0:
            raise FormatError("It don't look like a site bed6 file!")
        else:
            pass

    @classmethod
    def readbedfromregion(cls, file, left, right, label=None, genome=None, sample=None):
        """
        read region bed file and return a site like lst

        :param file:
        :param genome:
        :param sample:
        :return:
        """

        __res = []
        with open(file) as fh:
            for line in fh:
                chrom, s, e, _, _, strand = line.strip().split('\t')

                # skip the unusual chromosome
                if len(chrom) > 6:
                    continue

                if s == e:
                    continue
                _site = str(random.randrange(int(s), int(e)))

                __res.append('\t'.join([chrom,
                                        _site,
                                        _site,
                                        '.',
                                        '.',
                                        strand]))
        if sample:
            __res = random.sample(__res, sample)

        return cls(__res, left, right, label, genome)

    def write(self, file):
        self.mtx.to_csv(file, sep=',', index=False)

    @staticmethod
    def randomfasta(fasta, maxlength, sample, label, fileout=None):
        """
        random select sequence from the fasta file
        !Warning! ONLY SUPPORT THE STRANDARDS ENSEMBL mRNA FASTA file

        :param fasta:
        :param maxlength:
        :param sample:
        :return:
        """
        _id = re.compile(r'>(.*?) ')
        _gt = re.compile(r'gene_biotype:(.*?) ')
        _mid = maxlength // 2

        if fasta.endswith('gz'):
            fh = gzip.open(fasta, 'rb')
        else:
            fh = open(fasta)

        res = defaultdict(str)

        skip = False
        currentid = ''

        for line in fh:
            try:
                line = line.decode('utf-8')
            except:
                pass

            if line.startswith(">"):
                geneid = re.findall(_id, line)
                genetype = re.findall(_gt, line)

                if not geneid and not genetype:
                    raise FormatError("Fasta file don't look like strandard ensembl fasta file!")

                if genetype[0] != 'protein_coding':
                    skip = True
                else:
                    skip = False
                currentid = geneid[0]
            else:
                if not skip:
                    res[currentid] += line.strip()
        geneid = list(res.keys())

        _seq = []
        for i in range(sample):
            seq = res[geneid[random.choice(list(range(len(geneid))))]]

            # if len(seq) <= maxlength:
            #     continue
            try:
                index = random.randrange(_mid, len(seq) - _mid)
                seq = seq[index - _mid:index + _mid]
            except:
                continue

            if "N" in seq:
                continue
            _seq.append(seq)

        mtx = pd.DataFrame(_seq, columns=['sequence'])
        if label != None:
            mtx['label'] = label
        else:
            mtx['label'] = 'NaN'

        if fileout:
            mtx.to_csv(fileout, sep=',', index=False)
        else:
            return mtx

    @staticmethod
    def revseq(seq):

        seqdic = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "N": "N"
        }
        res = ''.join([seqdic[i] for i in list(seq)])
        return res


if __name__ == '__main__':
    a = SeqCrawler.readbed('splicesite.bed', 50, 50, 1, '/mnt/raid61/Microwell/mm10/fasta/genome.fa')
    a.write('test.csv')
