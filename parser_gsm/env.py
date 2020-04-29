import os, sys
import django
sys.path.append('/mnt/Storage2/home/zhengrongbin/project/geo_parser/py3/dc')
sys.path = sys.path[::-1]
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "dc.settings")
django.setup()
from datacollection import models
