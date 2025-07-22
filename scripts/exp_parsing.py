import pandas as pd
from pyparsing import (
    Word, alphanums, infixNotation, opAssoc, Keyword,
    ParserElement, ParseResults
)
from typing import Union, List, Any, Mapping

class BooleanExpressionParser:
    """
    A parser for KEGG boolean expressions using pyparsing.

    Handles expressions of the following form:
    - `gene1 AND gene2`
    - `gene1 OR gene2`
    - `NOT gene1`
    - Nested expressions like `(gene1 AND gene2) OR NOT gene3`
    """
    def __init__(self):
        # Enable performance optimization
        ParserElement.enablePackrat()

        # Define the parser
        self.GENE = Word(alphanums + "_:.-")
        self.AND = Keyword("AND")
        self.OR = Keyword("OR")
        self.NOT = Keyword("NOT")

        self.GENE_EXPRESSION = infixNotation(self.GENE,
            [
                (self.NOT, 1, opAssoc.RIGHT),
                (self.AND, 2, opAssoc.LEFT),
                (self.OR,  2, opAssoc.LEFT),
            ])

    def parse_expression(self, expression: str) -> ParseResults:
        """
        Parse a boolean expression string into a structured format.
        
        :param expression: The boolean expression to parse.
        :return: A ParseResults object representing the parsed expression.
        """
        return self.GENE_EXPRESSION.parseString(expression)

    def evaluate(self, parsed: Union[str, List[Any], ParseResults],
                 gene_row: Mapping[str, bool]) -> bool:
        """Evaluate a parsed boolean expression against a gene row.

        :param parsed: The parsed expression (from parse_expression).
        :param gene_row: A mapping of gene names to boolean values indicating presence.
        :return: True if the expression evaluates to true for the given gene row, False otherwise.
        """
        if isinstance(parsed, ParseResults):
            return self.evaluate(parsed.as_list()[0], gene_row)
        elif isinstance(parsed, str):
            return bool(gene_row.get(parsed, False))
        elif isinstance(parsed, list):
            if len(parsed) == 1:
                return self.evaluate(parsed[0], gene_row)
            elif parsed[0] == 'NOT':
                return not self.evaluate(parsed[1], gene_row)
            elif 'AND' in parsed:
                return all(self.evaluate(p, gene_row) for p in parsed if p != 'AND')
            elif 'OR' in parsed:
                return any(self.evaluate(p, gene_row) for p in parsed if p != 'OR')
        raise ValueError(f"Unexpected expression format: {parsed}")


# Unit testing
import unittest

class TestBooleanExpressions(unittest.TestCase):
    def setUp(self):
        self.parser = BooleanExpressionParser()

    def test_parse_expression(self):
        # Test parsing a simple expression
        parsed = self.parser.parse_expression("gene1 AND gene2 OR NOT gene3")
        self.assertIsInstance(parsed, ParseResults)
        flat_list = parsed.as_list(flatten=True)
        self.assertEqual(len(flat_list), 6)  # Should have 5 elements in the parsed result
        self.assertEqual(flat_list[0], 'gene1')
        self.assertEqual(flat_list[1], 'AND')
        self.assertEqual(flat_list[2], 'gene2')
        self.assertEqual(flat_list[3], 'OR')
        self.assertEqual(flat_list[4], 'NOT')
        self.assertEqual(flat_list[5], 'gene3')
    
    def test_evaluate_expression(self):
        # Test evaluating a simple expression
        gene_row = {
            'gene1': True,
            'gene2': False,
            'gene3': True,
            'gene4': False
        }

        strings_and_results = [
            ("gene1 AND gene2", False),
            ("gene1 OR gene2", True),
            ("NOT gene3", False),
            ("gene1 AND NOT gene2", True),
            ("gene3 OR (gene2 AND NOT gene1)", True),
            ("gene1 AND gene3 OR NOT gene2", True),
            ("(gene1 OR gene2) OR (gene3 AND NOT gene4)", True),
            ("(gene1 OR gene2) AND (gene3 AND gene4)", False),
            ("gene1 OR gene2 OR gene3 OR gene4", True),
            ("gene1 AND gene2 AND gene3 and gene4", False),
            ("(gene1 AND gene2) OR (gene3 AND gene4)", False),
            ("NOT (gene1 AND gene2)", True),
            ("NOT (gene1 OR gene2)", False),
        ]
        for expr_str, expected in strings_and_results:
            parsed = self.parser.parse_expression(expr_str)
            result = self.parser.evaluate(parsed, gene_row)
            self.assertEqual(result, expected, f"Failed for expression: {expr_str}")

if __name__ == "__main__":
    unittest.main()