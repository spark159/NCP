#!/usr/bin/env python

import sys
import math

class MarkovModel:
    def __init__(self, order):
        self.order = order

    # Train Markove Model from a list of DNA sequences
    def train(self, seq_list):
        # Return if no input
        if len(seq_list) <= 0:
            return

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

        # Normalize and take the log value of it
        for prob_order in self.prob:
            for prob_pos in prob_order:
                total_count = sum(prob_pos.values())
                total_count = float(total_count)
                for prev_cur, count in prob_pos.items():
                    prob_pos[prev_cur] = math.log(count / total_count)

    # Calculate log-probability of a given sequence
    def predict(self, seq):
        if self.seq_len != len(seq):
            return

        score = 0.0
        for i in range(self.seq_len):
            if i >= self.order:
                prev_cur = seq[i-self.order:i+1]
                if prev_cur not in self.prob[self.order][i]:
                    score += -100
                else:
                    score += self.prob[self.order][i][prev_cur]
            else:
                prev_cur = seq[:i+1]
                if prev_cur not in self.prob[i][i]:
                    score += -100
                else:
                    score += self.prob[i][i][prev_cur]
        return score

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
