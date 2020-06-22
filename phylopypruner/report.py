"""
Routines for formatting output.
"""

from __future__ import print_function
import sys
import os
from textwrap import wrap
# ensures that input is working across Python 2.7 and 3+
try:
    import __builtin__
    input = getattr(__builtin__, 'raw_input')
except (ImportError, AttributeError):
    pass


def supports_color():
    """
    Returns True if the running system's terminal supports color, and False
    otherwise.
    """
    plat = sys.platform
    supported_platform = plat != "Pocket PC" and (plat != "win32" or
                                                  "ANSICON" in os.environ)
    # isatty is not always implemented, #6223.
    is_a_tty = hasattr(sys.stdout, "isatty") and sys.stdout.isatty()
    if not supported_platform or not is_a_tty:
        return False
    return True


if supports_color():
    COLOR_SUPPORT = True
else:
    COLOR_SUPPORT = False
NORMAL = "\033[0m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
PURPLE = "\033[35m"
CYAN = "\033[36m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"


def warning(message, display=True):
    """Write the provided message, with the text 'warning: ' in front of it, to
    stderr.

    Parameters
    ----------
    message : str
        This is the string you want to print.

    display : bool
        Write this message to stderr, if display is True.

    Returns
    -------
    warning_message : str
        The formatted string which is also printed to stderr.
    """
    warning_message = "{}warning: {}{}".format(YELLOW, NORMAL, message)
    if display:
        print(warning_message, file=sys.stderr)
    return warning_message


def error(message, display=True):
    """Write the provided message, with the text 'error: ', in front of it to
    stderr.

    Parameters
    ----------
    message : str
        This is the string you want to print.

    display : bool
        Write this message to stderr, if display is True.

    Returns
    -------
    error_message : str
        The formatted string which is also printed to stderr.
    """
    error_message = "{}error: {}{}".format(RED, NORMAL, message)
    if display:
        print(error_message, file=sys.stderr)
    return error_message


def tip(message, display=True):
    """Write the provided message, with the text 'tip: ' in front of it, to
    stderr.

    Parameters
    ----------
    message : str
        This is the string you want to print.

    display : bool
        Write this message to stderr, if display is True.

    Returns
    -------
    warning_message : str
        The formatted string which is also printed to stderr.
    """
    tip_message = "{}tip: {}{}".format(GREEN, NORMAL, message)
    if display:
        print(tip_message, file=sys.stderr)
    return tip_message


def yes_or_no(question):
    """Takes a question as an input and prompts the user for a yes or no. Returns
    True if the answer is yes and False if the answer is no.

    Parameters
    ----------
    question : str
        This text is displayed, before asking for a yes or a no.

    Returns
    -------
    True or False
        True if answer is 'y', False if answer is 'n'.
    """
    # make input work the same way in both Python 2 and 3
    answer = input(question + " (y/n): ".lower().rstrip())
    while not answer in ("y", "n"):
        answer = input("please answer yes or no" + " (y/n): ".lower().rstrip())
    return answer[0] == "y"


def progress_bar(message, replace=True, display=True):
    """Takes a text string as an input and prints it to standard error, with a
    greater than sign '>' prepended to the text. Only returns the message
    itself, if the boolean display is set to False (DEFAULT: True). By default,
    the last line is replaced as this function is intended to be used to
    display a continues progress. If you do not wish to overwrite the last
    line, then set replace to False.

    Parameters
    ----------
    message : str
        This is the string you want to print.

    replace : bool
        Overwrite the last line as in a progress bar.

    display : bool
        Write this message to stderr, if display is True.

    Returns
    -------
    progress : str
        Your message with the string '->' prepended to it.
    """
    progress = "{} {}".format(colorize(u">", "green"), message)
    if display and replace:
        sys.stdout.flush()
        print(progress, file=sys.stderr, end="\r")
    elif display:
        print(progress, file=sys.stderr)
    return progress


def print_path(path, display=True):
    """Takes the path to a file as an input and prints the given path in blue.

    Parameters
    ----------
    path : str
      The path you wish to colorize.

    display : bool
      Print the results to stderr if this variable is True. Defaults to True.

    Returns
    -------
    path_in_blue : str
      The provided path in blue.
    """
    if display:
        print(path, file=sys.stderr)
    return path


def underline(text):
    """Returns the provided text with an underline.

    Parameters
    ----------
    text : str
        The text you wish to underline.

    Returns
    -------
    underlined_text : str
        The provided text, underlined.
    """
    underlined_text = ("{}{}{}".format(UNDERLINE, text, NORMAL))
    return underlined_text


def bold(text):
    """Returns the provided text in bold.

    Parameters
    ----------
    text : str
        The text you wish to make bold.

    Returns
    -------
    bold_text : str
        The provided text in bold.
    """
    bold_text = ("{}{}{}".format(BOLD, text, NORMAL))
    return bold_text


def display_otus(otus, display=True):
    """Takes a list of OTUs as an input and returns the OTUs as a string, where
    each OTU, except for the last OTU, is separated by a comma. OTUs are also
    displayed in the color red.

    Parameters
    ----------
    otus : list
      A list of OTUs.

    display : bool
      The OTUs are printed to stderr if display is set to True. True by
      default.

    Returns
    -------
    formatted_otus : str
      A list of the OTUs as a text string.
    """
    formatted_otus = ", ".join(
        [colorize(otu, "red") for otu in otus])
    formatted_otus = "\n    ".join(wrap(formatted_otus, 150))
    if display:
        print("    " + formatted_otus, file=sys.stderr)
    return "    " + formatted_otus


def colorize(text, color):
    """Returns the provided text in the chosen color. The color you pick should
    be one of the following: "red", "green", "yellow", "blue", "purle", or
    "cyan". If colors are not supported by your terminal, then no color is
    applied.

    Parameters
    ----------
    text : str
      The text you wish to colorize.

    color : str
      The name of the color you wish to use, in lowercase.

    Returns
    -------
    colored_text : str
      The provided text in color, given that colors are supported.
    """
    palette = {"red": RED,
               "green": GREEN,
               "yellow": YELLOW,
               "blue": BLUE,
               "purple": PURPLE,
               "cyan": CYAN}
    colored_text = "{}{}{}".format(palette[color], text, NORMAL)
    if COLOR_SUPPORT:
        return colored_text
    return text


def format_otus(otus):
    """Returns the provided OTUs as a formatted string, where each OTU in the
    list is separated by a comma and a space (', '). Also add a newline and 8
    spaces to string. If the list of OTUs is empty, then return 'None'.

    Parameters
    ----------
    otus : list
      List of the OTUs.

    Returns
    -------
    otus_as_str : str
      Formatted string of the provided OTUs or 'None' if OTUs is empty.
    """
    if not otus:
        otus_as_str = "None"
    else:
        otus_as_str = "\n" + " " * 8 + ", ".join([otu for otu in otus])
    return otus_as_str
