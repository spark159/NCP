#!/usr/bin/env python

import sys
import math
import numpy as np

class MarkovModel:
    def __init__(self, order):
        self.order = order

    # Train Markove Model from genome and a list of DNA sequences
    def train(self, genome, seq_list):
        # Return if no input
        if len(seq_list) <= 0:
            return

        # Compute background signal
        self.bg_prob = [{} for i in range(self.order+1)]
        for seq in genome.values():
            for i in range(len(seq)):
                for j in range(self.order + 1):
                    prev = seq[i-j:i]
                    cur = seq[i]
                    assert cur in "ACGT"
                    prev_cur = prev + cur
                    prob_j = self.bg_prob[j]
                    if prev_cur not in prob_j:
                        prob_j[prev_cur] = 1
                    else:
                        prob_j[prev_cur] += 1

        # Normalize and take the log values of them
        for prob_order in self.bg_prob:
            total_count = sum(prob_order.values())
            total_count = float(total_count)
            for prev_cur, count in prob_order.items():
                prob_order[prev_cur] = math.log(count / total_count)

        # Make sure that sequences are of the same length
        self.seq_len = len(seq_list[0])
        for seq in seq_list[1:]:
            assert self.seq_len == len(seq)

        # self.prob[order][position][previous_bases + current_base]
        self.prob = [[{} for j in range(self.seq_len)] for i in range(self.order+1)]
        for seq in seq_list:
            for i in range(self.seq_len):
                for j in range(self.order + 1):
                    prev = seq[i-j:i]
                    cur = seq[i]
                    assert cur in "ACGT"
                    prev_cur = prev + cur
                    prob_i = self.prob[j][i]
                    if prev_cur not in prob_i:
                        prob_i[prev_cur] = 1
                    else:
                        prob_i[prev_cur] += 1

        # Normalize and take the log values of them
        for prob_order in self.prob:
            for prob_pos in prob_order:
                total_count = sum(prob_pos.values())
                total_count = float(total_count)
                for prev_cur, count in prob_pos.items():
                    prob_pos[prev_cur] = math.log(count / total_count)

    # Calculate log-probability of a given sequence
    def predict(self, seq, background = True):
        if self.seq_len != len(seq):
            return

        score = 0.0
        for i in range(self.seq_len):
            if i >= self.order:
                prev_cur = seq[i-self.order:i+1]
                order = self.order
            else:
                prev_cur = seq[:i+1]
                order = i

            if prev_cur not in self.prob[order][i]:
                score = (-sys.float_info.max / 2)
            else:
                score += self.prob[order][i][prev_cur]

            if background:
                if prev_cur not in self.bg_prob[order]:
                    score = (-sys.float_info.max / 2)
                else:
                    score += -self.bg_prob[order][prev_cur]
                
        return score

    # Calculate the sum of weight-factor of given sub sequence (forward)
    def get_forward(self, seq):
        N = len(seq)
        F = [1.0] * (N + 1)
        for i in range(N):
            bp = seq[i]
            bg_prob = math.exp(self.bg_prob[0][bp])
            F[i+1] = F[i] * bg_prob
            if i >= self.seq_len:
                F[i+1] += F[i+1-self.seq_len] * math.exp(self.predict(seq[i-self.seq_len+1:i+1], False))
        return F

    # Calculate the sum of weight-factor of given sub sequence (reverse)
    def get_reverse(self, seq):
        N = len(seq)
        R = [1.0] * (N + 1)
        for i in reversed(range(N)):
            bp = seq[i]
            bg_prob = math.exp(self.bg_prob[0][bp])
            R[i] = R[i+1] * bg_prob
            if i + self.seq_len < N:
                R[i] += R[i + self.seq_len] * math.exp(self.predict(seq[i:i+self.seq_len], False))
        return R

    # calculate the probabilty of NCP would be on ith position of given sequence
    def prob_i(self, seq, i, F, R):
        N = len(seq)
        i -= (self.seq_len / 2)
        if i < 0 or i + self.seq_len > N:
            return 0.0
        return F[i] * math.exp(self.predict(seq[i:i+self.seq_len], False)) * R[i + self.seq_len] / R[0]

    # Display some detail about the model
    def help(self):
        print >> sys.stderr, "%d-order Markov Model" % (self.order)
        for order in range(self.order + 1):
            for pos in range(len(self.prob[order])):
                if pos % 20 != 0:
                    continue
                if len(self.prob[order][pos]) == 0:
                    continue
                print >> sys.stderr, "%d-order at pos %d" % (order, pos)
                print >> sys.stderr, "\t", self.prob[order][pos]


# To be implemented
class InterpolatedMarkovModel:
    def __init__(self):
        None
        

# To be implemented
class HiddenMarkovModel:
    def __init__(self):
        None
