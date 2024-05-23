
import numpy as np
from typing import Literal


class GaussProt(object):

    def __init__(
            self,
            schema: dict,
            standardize_schema=True,
            shard_size: int = 10000,
            bandwidth: float = 0.8,
            model_type: Literal['discrete', 'continuous'] = 'discrete',
            discrete_model_length: int = None,
            padded: bool = True,
            verbose: bool = False
    ):
        """
        GaussProt main class.


        :param schema:
        :param standardize_schema:
        :param shard_size:
        :param bandwidth:
        :param model_type:
        :param discrete_model_length:
        :param padded:
        :param verbose:
        """
        self.standardize_schema = standardize_schema
        self.shard_size = shard_size
        self.bandwidth = bandwidth
        self.model_type = model_type
        self.verbose = verbose
        if self.model_type not in ['discrete', 'continuous']:
            raise ValueError('model_type must be discrete or continuous')
        self.discrete_model_length = discrete_model_length
        self.padded = padded
        self.schema = schema
        self._validate_schema()
        if standardize_schema:
            self._standardize_schema()

    def generate_models(self, sequences: list[str]) -> np.ndarray:
        """
        Generate models of sequences.
        :param sequences: List of protein sequences
        :return: 2d array
        """
        if self.verbose:
            print(f'Generating models with parameters:')
            print(f'model_type: {self.model_type}')
            print(f'bandwidth: {self.bandwidth}')
            print(f'standardize_schema: {self.standardize_schema}')
            print(f'shard_size: {self.shard_size}')
            if self.model_type == 'discrete' and self.discrete_model_length:
                print(f'discrete_model_length: {self.discrete_model_length}')
            elif self.model_type == 'continuous':
                print(f'padded: {self.padded}')

        if len(sequences) <= self.shard_size:
            return self._generate_models(sequences)

        shards = []
        for i in range(0, len(sequences), self.shard_size):
            slice_ = sequences[i: i+self.shard_size]
            if len(slice_) < 1:
                break
            shards.append(self._generate_models(slice_))
        return np.concatenate(shards)

    def _generate_models(self, sequences: list[str]) -> np.ndarray:
        sequence_lengths = [len(seq) for seq in sequences]
        ML = max(sequence_lengths)
        N = len(sequences)
        vector_length = self.discrete_model_length or 10
        MVL = vector_length * (ML + 2)
        weights_frame = np.zeros((N, ML))
        matrix_frame = np.zeros((N, MVL))
        x = np.linspace(-2.3263, 2.3263, 3 * vector_length)
        pdf_template = self._pdf(x)
        # vectorize weights
        for i, seq in enumerate(sequences):
            weights = [self.schema[a] for a in seq]
            weights_frame[i, 0:len(weights)] = weights
        # compute continuous profiles
        idx = list(zip(range(ML), range(0, MVL - vector_length, vector_length)))
        for wi, si in idx:
            ith_frame = np.zeros_like(matrix_frame)
            ith_frame[:, si:si + 3 * vector_length] = pdf_template[np.newaxis, :] * weights_frame[:, wi][:, np.newaxis]
            matrix_frame += ith_frame
        if self.model_type == 'continuous':
            if self.padded:
                return matrix_frame
            non_padded = np.array([])
            for i in range(N):
                SL = sequence_lengths[i]
                CVL = vector_length * (SL + 2)
                slice = matrix_frame[i, vector_length: CVL - vector_length]
                non_padded[i] = slice
            return non_padded
        # compute discrete profiles
        if self.discrete_model_length:
            discrete_matrix_frame = np.zeros((N, self.discrete_model_length))
        else:
            discrete_matrix_frame = np.array([])
        for i in range(N):
            SL = sequence_lengths[i]
            CVL = vector_length * (SL + 2)
            slice = matrix_frame[i, vector_length: CVL - vector_length]
            if self.discrete_model_length:
                windows = np.split(slice, vector_length)
                discrete_matrix_frame[i, :] = self._auc(windows)
            else:
                windows = np.split(slice, SL)
                discrete_matrix_frame[i] = self._auc(windows)
        return discrete_matrix_frame

    def _pdf(self, x):
        y = np.exp(-x ** 2 / (2 * self.bandwidth)) / (self.bandwidth * np.sqrt(2 * np.pi))
        return (y - y.min()) / (y.max() - y.min())

    @staticmethod
    def _auc(windows):
        return np.array([array.sum() / 2.0 for array in windows])

    def _validate_schema(self):
        valid_letters = [
            'C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M'
        ]
        if not all([k in self.schema for k in valid_letters]):
            raise ValueError('Invalid schema.\nValid letters: {}'.format(valid_letters))

    def _standardize_schema(self):
        standard_params = {}
        vals = list(self.schema.values())
        min_ = min(vals)
        max_ = max(vals)
        for k, v in self.schema.items():
            self.schema[k] = 2 * (v - min_) / (max_ - min_) - 1

