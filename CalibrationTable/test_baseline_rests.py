import unittest
from unittest.mock import patch

import resume_baseline


class BaselineRestTests(unittest.TestCase):
    def test_only_scheduled_rests_change(self):
        original = resume_baseline.baseline.core.Engine
        shortened = resume_baseline.engine_with_rests(5, 30)
        engine = object.__new__(shortened)
        with patch.object(original, 'wait') as wait:
            engine.wait(30, 'scheduled rest after new synthesis')
            wait.assert_called_with(5, 'scheduled rest after new synthesis')
            engine.wait(300, 'scheduled rest after new synthesis')
            wait.assert_called_with(30, 'scheduled rest after new synthesis')
            engine.wait(60, 'waiting for system memory')
            wait.assert_called_with(60, 'waiting for system memory')
        self.assertIs(shortened.solve, original.solve)
        self.assertIs(shortened.compute, original.compute)


if __name__ == '__main__':
    unittest.main()
