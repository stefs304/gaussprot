
import unittest
import numpy as np
from gaussprot import GaussProt
from .data import (
    VALID_SCHEMA,
    INVALID_SCHEMA,
    VALID_SEQUENCES,
    INVALID_SEQUENCES
)


class TestGaussProt(unittest.TestCase):

    def test_discrete_case(self):
        gp = GaussProt(schema=VALID_SCHEMA, model_type='discrete')
        models = gp.generate_models(VALID_SEQUENCES)
        self.assertEquals(len(models), len(VALID_SEQUENCES))
        self.assertEquals(type(models), list)
        for i in range(len(VALID_SEQUENCES)):
            self.assertEquals(len(VALID_SEQUENCES[i]), len(models[i]))

        fixed_length = 5
        gp = GaussProt(schema=VALID_SCHEMA, model_type='discrete', discrete_model_length=fixed_length)
        models = gp.generate_models(VALID_SEQUENCES)
        self.assertEquals(len(models), len(VALID_SEQUENCES))
        self.assertEquals(type(models), list)
        for i in range(len(VALID_SEQUENCES)):
            self.assertEquals(fixed_length, len(models[i]))

    def test_continuous_case(self):
        gp = GaussProt(schema=VALID_SCHEMA, model_type='continuous')
        models = gp.generate_models(VALID_SEQUENCES)
        self.assertEquals(len(models), len(VALID_SEQUENCES))
        self.assertEquals(type(models), list)
        for i in range(len(VALID_SEQUENCES)):
            expected_length = max([len(seq) for seq in VALID_SEQUENCES]) * gp.multiplier
            self.assertEquals(expected_length, len(models[i]))

        gp = GaussProt(schema=VALID_SCHEMA, model_type='continuous', padded=False)
        models = gp.generate_models(VALID_SEQUENCES)
        self.assertEquals(len(models), len(VALID_SEQUENCES))
        self.assertEquals(type(models), list)
        for i in range(len(VALID_SEQUENCES)):
            expected_length = len(VALID_SEQUENCES[i]) * gp.multiplier
            self.assertEquals(expected_length, len(models[i]))

    def test_exceptions(self):

        with self.assertRaises(ValueError):
            GaussProt(schema=INVALID_SCHEMA)

        gp = GaussProt(schema=VALID_SCHEMA)
        with self.assertRaises(KeyError):
            gp.generate_models(sequences=INVALID_SEQUENCES)

    def test_init_signature(self):
        raised = False
        try:
            GaussProt(
                schema=VALID_SCHEMA,
                standardize_schema=True,
                shard_size=1000,
                bandwidth=0.5,
                model_type='discrete',
                discrete_model_length=11,
                verbose=True,
                validate_schema=True
            )
            GaussProt(
                schema=VALID_SCHEMA,
                standardize_schema=False,
                shard_size=1000,
                bandwidth=0.1,
                model_type='continuous',
                padded=False,
                verbose=False,
                validate_schema=False
            )
        except Exception as e:
            print(e)
            raised = True
        self.assertFalse(raised, 'v <= 0.0.1 signature test failed.')

    def test_simply_encode(self):
        gp = GaussProt(schema=VALID_SCHEMA, model_type='discrete', padded=False)
        encoded = gp.simply_encode(VALID_SEQUENCES)
        self.assertEqual(len(encoded), len(VALID_SEQUENCES))
        self.assertEqual(type(encoded), list)
        for i in range(len(VALID_SEQUENCES)):
            self.assertEqual(len(VALID_SEQUENCES[i]), len(encoded[i]))

        gp = GaussProt(schema=VALID_SCHEMA, model_type='discrete', padded=True)
        encoded = gp.simply_encode(VALID_SEQUENCES)
        max_len = max([len(seq) for seq in VALID_SEQUENCES])
        self.assertEqual(len(VALID_SEQUENCES), len(encoded))
        self.assertEqual(type(encoded), list)
        for enc_seq in encoded:
            self.assertEqual(max_len, len(enc_seq))

