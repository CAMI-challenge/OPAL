Feature: Run Scripts without error

  Background: Initial Setup

  Scenario: Run OPAL with basic command
    Given I create the directory "input"
    And I create the directory "output"
    And I copy the example data files
      | source                                     | dest        |
      | goldstandard_low_1.bin     | input |
      | test1.profile      | input |
      | test2.profile      | input |
    When I run the command
    """
    ../opal.py -g input/goldstandard_low_1.bin -o output input/test1.profile input/test2.profile
    """
    Then the exit code should be 0
