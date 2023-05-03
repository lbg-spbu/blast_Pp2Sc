import re
from functools import wraps
from time import sleep
from typing import Optional


def retries(times: int = 3, delay: int = 0, exception=Exception):
    # TODO description
    """

    """

    def decorator_function(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            exception_message = ""

            for _ in range(times):
                try:
                    result = func(*args, **kwargs)
                except exception as e:
                    sleep(delay)
                    exception_message = str(e)
                    continue
                else:
                    return result

            raise exception(exception_message)

        return wrapper

    return decorator_function


def parse_gene_name(description: str) -> Optional[str]:
    pattern = r"(?P<gene_name>[a-zA-Z]{3}\d{1,2})"
    try:
        res = re.search(pattern=pattern, string=description)
    except Exception:
        return None

    if res is None:
        return None

    return res.group()
