def is_string_like(obj):
    try:
        obj + ''
    except:
        return False
    return True
