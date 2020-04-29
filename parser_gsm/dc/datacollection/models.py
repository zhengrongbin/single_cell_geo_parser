from django.db import models
from django.contrib.auth.models import User
from django.utils.encoding import smart_str
import os
import django.db.models.options as options
options.DEFAULT_NAMES = options.DEFAULT_NAMES + ('in_db',)
from django.db.models.signals import m2m_changed


try:
    import json
except ImportError:
    import simplejson as json

#drop?
SPECIES_CHOICES = (
    ('hg', 'homo sapien'),
    ('mm', 'mus musculus'),
)
#drop?
EXPERIMENT_TYPE_CHOICES = (
    ('chip', 'ChIP-Chip'),
    ('seq', 'ChIP-Seq'),
)

SUBMISSION_STATUS = (
    ('pending', 'Pending'),
    ('closed', 'Imported/Closed'),
    ('n/a', 'Not Appropriate'),
)

PAPER_STATUS = (
    ('imported', 'paper entered awaiting datasets'),
    ('datasets', 'datasets imported awaiting download'),
    ('transfer', 'datasets download in progress'),
    ('downloaded', 'datasets downloaded awaiting analysis'),
    ('complete', 'analysis complete/complete'),
    ('error', 'error/hold- see comments'),
)

GENERAL_STATUS = (
    ('new', 'newly imported'),
    ('ok', 'validated'),
)

DATASET_STATUS = (
    ('inherited', 'inherited from DC 1.0'),
    ('validated', 'meta-info validated awaiting file download'),
    ('auto-parsed', 'meta-info extracted automatically awaiting validation'),
    ('transfer', 'file is downloading'),
    ('downloaded', 'downloaded/closed'),
    ('complete', 'analysis complete'),
    ('error', 'error/hold- see comments'),
    ('?', 'not sure')
)

SAMPLE_STATUS = (
    ('new', 'sample created'),
    ('checked', 'sample checked, awaiting importation'),
    ('inherited', 'sample inherited from DC 1.0'),
    ('running', 'analysis is running, awaiting completion'),
    ('complete', 'analysis complete'),
    ('ignored', 'ignored or wrong sample'),
)

TEAMS = (
    ('admin', 'Administrators'),
    ('paper', 'paper collection team'),
    ('data', 'data collection team'),
)

FACTOR_TYPES = (
    ('tf', 'transcription factor'),
    ('hm', 'histone mark'),
    ('cr', 'chromatin regulator'),
    ('ca','chromatin accessibility'),
    ('predicted transcription factor', 'predicted transcription factor'),
    ('predicted chromatin regulator', 'predicted chromatin regulator'),
    ('both predicted transcription factor and chromatin regulator',
     'both predicted transcription factor and chromatin regulator'),
    ('not sure','not sure'),
    ('other', 'other'),
)

#pending is the default submission status
DEFAULT_SUBMISSION_STATUS = SUBMISSION_STATUS[0][0]


class DCModel(models.Model):
    """Implements common fns that will be inherited by all of the models
    NOTE: this is an ABSTRACT class, otherwise django will try to make a db
    for it.
    """

    class Meta:
        abstract = True

    def to_json(self):
        """Returns the model as a json string"""
        tmp = {}
        for k in list(self.__dict__.keys()):
            tmp[k] = "%s" % self.__dict__[k]
        return json.dumps(tmp)


class Papers(DCModel):
    """Papers are publications that are publicly available
    The fields are:
    pmid - the pubmed id of the publication
    unique_id - the unique identifier for an external db, e.g. GSEID
    user - the user who currated the dataset
    title - the title of the publication
    reference - shorthand of the publication details
    abstract - the abstract of the publication
    pub_date - the publication date
    authors - list of the paper's authors
    last_aut_email - the email address of the last/corresponding author
    journal- paper's journal

    status- SEE PAPER_STATUS above
    comments- any comments about this paper
    """
    #def __init__(self, *args):
    #    super(Papers, self).__init__(*args)
    #    self._meta._donotSerialize = ['user']
    class Meta:
        verbose_name_plural = 'papers'

    pmid = models.IntegerField(null=True, blank=True, default=None)
    #NOTE: papers can have multiple unique_ids attached--if so, comma-sep them
    unique_id = models.CharField(max_length=255, null=True, blank=True, default="")
    user = models.ForeignKey(User, null=True, blank=True, default=None, on_delete=models.CASCADE)
    title = models.CharField(max_length=255, null=True, blank=True, default="")
    reference = models.CharField(max_length=255, null=True, blank=True, default="")
    abstract = models.TextField(null=True, blank=True, default="")
    pub_date = models.DateField(null=True, blank=True, default=None)
    date_collected = models.DateTimeField(null=True, blank=True, default=None)
    authors = models.CharField(max_length=1000, null=True, blank=True, default="")
    last_auth_email = models.EmailField(null=True, blank=True, default=None)

    journal = models.ForeignKey('Journals', on_delete=models.CASCADE,
                                null=True, blank=True, default=None)
    status = models.CharField(max_length=255, choices=PAPER_STATUS,
                              null=True, blank=True, default="imported")
    #a place for curators to add comments

    comments = models.TextField(null=True, blank=True, default="")

    def _sampleAggregator(sample_field):
        """Given a sample field, tries to aggregate all of the associated
        samples"""
        #pluralizing the dset_field
        #exceptions first
        if sample_field == 'spe':
            plural = 'species'
        elif sample_field == 'assembly':
            plural = 'assemblies'
        else:
            plural = sample_field + 's'

        def nameless(self):
            "Returns a list of %s associates with the papers" % plural
            ids = []
            vals = []

            samples = Samples.objects.filter(paper=self.id)
            for s in samples:
                val = getattr(s, sample_field)
                if val and val.id not in ids:
                    ids.append(val.id)
                    vals.append(smart_str(val))
            return vals

        return nameless

    #NOTE: these are killing the initial cache!!! AND we're only using one
    # of these--species.  For efficiency sake, commenting the rest out!
    #NOTE: factors, cell_lines, cts, cps, tts is also needed for search
    factors = property(_sampleAggregator('factor'))
    # platforms = property(_sampleAggregator('platform'))
    species = property(_sampleAggregator('species'))
    # assemblies = property(_sampleAggregator('assembly'))
    cell_types = property(_sampleAggregator('cell_type'))
    cell_lines = property(_sampleAggregator('cell_line'))
    cell_pops = property(_sampleAggregator('cell_pop'))
    tissue_types = property(_sampleAggregator('tissue_type'))
    # strains = property(_sampleAggregator('strain'))
    # conditions = property(_sampleAggregator('condition'))
    # disease_states = property(_sampleAggregator('disease_state'))

    def _aggUniqueIds(self):
        """aggregates the unique ids of the samples associated w/ the paper"""
        samples = Samples.objects.filter(paper=self.id)
        #NOTE: we are wrapping up the unique_id and the associated factor!
        return [[smart_str(s.unique_id), smart_str(s.factor)] for s in samples]

    sample_unique_ids = property(_aggUniqueIds)

    def _get_lab(self):
        """Returns the last author in the authors list"""
        try:
            return smart_str(self.authors.split(",")[-1:][0]).strip()
        except:
            return smart_str(self.authors).strip()


    lab = models.CharField(max_length=1000, null=True, blank=True, default="")

    pub_summary = models.CharField(max_length=1000, null=True, blank=True, default="")
    # a dirty trick to print first author, journal, and pub_date

    def __str__(self):
        return smart_str(self.title)

    def save(self):
        first_author = self.authors.split(",")[0] + ", et al."
        title = self.title
        journal = self.journal.name if self.journal else ""
        year = str(self.pub_date.year) if self.pub_date else ""
        self.reference = " ".join([first_author , title , journal, year])
        self.lab = self._get_lab()
        self.pub_summary = " ".join([first_author , journal, year])
        super(Papers, self).save()

#NOTE: _donotSerialize fields are not enumerated as records, just as keys
Papers._meta._donotSerialize = ['user']

#Dataset fields which we will aggregate and make into virtual paper fields
#NOTE: the only ones used by jscript is species and sample_unique_ids!!
Papers._meta._virtualfields = [#'lab', 'factors', 'platforms', 'species',
                               #'assemblies', 'cell_types', 'cell_lines',
                               #'cell_pops', 'tissue_types', 'strains',
                               #'conditions', 'disease_states',
                               'species',
                               'sample_unique_ids',
                               'factors', 'cell_types', 'cell_lines',
                               'cell_pops', 'tissue_types']




class Datasets(DCModel):
    """Datasets are SETS of samples, typically a set (1 or more) of treatment
    samples and a set (0 or more) of control samples.

    **They are the HIGHER-level representations; samples are the lower-level.

    Dataset fields:
    user_id - the user who currated the dataset
    paper_id - the paper that the dataset is associated with

    [Virtual fields- meta data inferred from treatment (and sometimes
    control) samples]
    [File fields]
    peak_bed
    peak_wig
    ...etc..
    """

    class Meta:
        verbose_name_plural = 'datasets'

    #MIG_NOTE: CAN'T depend on GSMIDS
    #def upload_factory(sub_dir):
    #    """a factory for generating upload_to_path fns--e.g. use to generate
    #    the various sub-directories we use to store the info associated w/
    #    a sample"""

    #    def upload_to_path(self, filename):
    #        """Returns the upload_to path for this dataset.
    #        NOTE: we are going to store the files by gsmid, e.g. GSM566957
    #        is going to be stored in: data/datasets/gsm566/957.
    #        I'm not sure if this is the place to validate gsmids, but it maybe
    #        """
    #        #20120821 THIS WILL NOT WORK: gsmid DNE!
    #        return os.path.join('data', 'datasets', 'gsm%s' % self.gsmid[3:6],
    #                            self.gsmid[6:], sub_dir, filename)

    #    return upload_to_path

    #def __init__(self, *args):
    #    super(Datasets, self).__init__(*args)
    #    self._meta._donotSerialize = ['user']

    @staticmethod
    def sample_filter(**params):
        """Returns a list of datasets whose associated samples have that
        information.
        NOTE: paper is a dataset field and a sample field--
              sample_filter will go to the paper level (for consistency)
              IF you want the datasets assoc. with a paper, use
              Datasets.objects.filter(paper = )
        """
        ret = []
        dsets = Datasets.objects.all()
        for d in dsets:
            valid = True
            for k in list(params.keys()):
                if getattr(d, k) != params[k]:
                    valid = False
                    break;
            if valid:
                ret.append(d)
        return ret

    def _sampleAggregator(sample_field):
        """given a sample field, like 'Factor'--TRIES:
        1. to take the factor of the first treatment if that is valid
        2. take the factor of the first control
        3. Otherwise None
        """

        def nameless(self):
            """Assumes that the assoc. sample %s entries are all consistent
            so this fn returns the first valid %s--otherwise None""" % \
            (sample_field, sample_field)
            val = None
            if self.treats:

                all_treat = self.treats.all()
                if all_treat:
                    val = getattr(all_treat[0], sample_field)
            elif self.conts:
                all_cont = self.conts.all()
                if all_cont:
                    val = getattr(all_cont[0], sample_field)
            if val:
                return val.name
            else:
                return val

        return nameless

    #user = the person curated/created the dataset
    user = models.ForeignKey(User, null=True, blank=True, default=None, on_delete=models.CASCADE)
    paper = models.ForeignKey('Papers', null=True, blank=True, default=None, on_delete=models.SET_NULL)

    treats = models.ManyToManyField('Samples', related_name="TREATS")
    conts = models.ManyToManyField('Samples', related_name="CONTS", blank=True, default=None) #null=True, 

    #treatments = models.CommaSeparatedIntegerField(max_length=255, null=True,
    #                                               blank=True)
    #controls = models.CommaSeparatedIntegerField(max_length=255, null=True,
    #                                             blank=True)

    paper = models.ForeignKey('Papers', null=True, blank=True, default=None, on_delete=models.SET_NULL)
    factor = models.ForeignKey('Factors', null=True, blank=True, default=None,on_delete=models.SET_NULL)
    species = models.ForeignKey('Species',
                                null=True, blank=True, default=None,on_delete=models.SET_NULL)
    cell_type = models.ForeignKey('CellTypes',
                                  null=True, blank=True, default=None,on_delete=models.SET_NULL)
    cell_line = models.ForeignKey('CellLines',
                                  null=True, blank=True, default=None, on_delete=models.SET_NULL)
    cell_pop = models.ForeignKey('CellPops',
                                 null=True, blank=True, default=None,on_delete=models.SET_NULL)
    strain = models.ForeignKey('Strains',
                               null=True, blank=True, default=None,on_delete=models.SET_NULL)
    disease_state = models.ForeignKey('DiseaseStates',
                                      null=True, blank=True, default=None,on_delete=models.SET_NULL)
    tissue_type = models.ForeignKey('TissueTypes', null=True, blank=True,
                                    default=None,on_delete=models.SET_NULL)

    full_text = models.TextField(null=True, blank=True, default="")
    result_folder = models.CharField(max_length=255, null=True, blank=True, default="")

    #FileFields
    # peak_file = models.FileField(upload_to=upload_factory("peak"),
    #                              null=True, blank=True, max_length=1024)
    # peak_xls_file = models.FileField(upload_to=upload_factory("peak"),
    #                                  null=True, blank=True, max_length=1024)
    # summit_file = models.FileField(upload_to=upload_factory("summit"),
    #                               null=True, blank=True, max_length=1024)

    # treat_bw_file = models.FileField(upload_to=upload_factory("bw"),
    #                                 null=True, blank=True, max_length=1024)
    # cont_bw_file = models.FileField(upload_to=upload_factory("bw"),
    #                                 null=True, blank=True, max_length=1024)
    #
    # conservation_file = models.FileField(upload_to=upload_factory("conservation"),
    #                                      null=True, blank=True, max_length=1024)
    # conservation_r_file = models.FileField(upload_to=upload_factory("conservation"),
    #                                        null=True, blank=True, max_length=1024)
    #ceas_file = models.FileField(upload_to=upload_factory("ceas"),
    #                             null=True, blank=True, max_length=1024)
    #summary_file = models.FileField(upload_to=upload_factory("meta"),
    #                                null=True, blank=True, max_length=1024)
    #
    # #NOTE: even though dhs stats are saved in a table, we're going to store it
    # #in meta
    # dhs_file = models.FileField(upload_to=upload_factory("meta"),
    #                             null=True, blank=True, max_length=1024)

    date_created = models.DateTimeField(blank=True, null=True, default=None)

    status = models.CharField(max_length=255, choices=DATASET_STATUS,
                              default="new")
    comments = models.TextField(blank=True, default="")

    fastqc = models.CharField(max_length=255, null=True)
    mapped = models.CharField(max_length=255, null=True)
    map_ratio = models.CharField(max_length=255, null=True)
    pbc = models.CharField(max_length=255, null=True)
    peaks = models.CharField(max_length=255, null=True)
    frip = models.CharField(max_length=255, null=True)
    dhs = models.CharField(max_length=255, null=True)
    motif = models.CharField(max_length=255, null=True)
    qc_summary = models.CharField(max_length=255, null=True)
    def _printInfo(self):
        """Tries to print the treatment and controls list like this:
        GSMXXX,GSMYYY::GSMZZZ where the :: is a separator btwn the two lists"""
        treat = []
        control = []
        if self.treatments:
            treat = [Datasets.objects.get(id=d) \
                     for d in self.treatments.split(',')]
        if self.controls:
            controls = [Datasets.objects.get(id=d) \
                        for d in self.controls.split(',')]
        return "%s::%s" % (",".join(treat), ",".join(control))

    # For admin interface

    def journal_name(self):
        if self.paper:
            return self.paper.journal.name if self.paper.journal else None
        else:
            return None

    journal_name.admin_order_field = 'paper__journal__name'

    def journal_impact_factor(self):
        if self.paper:
            return self.paper.journal.impact_factor if self.paper.journal else None
        else:
            return None

    def paper_pmid(self):
        if self.paper:
            return self.paper.pmid
        else:
            return None

    def __str__(self):
        return smart_str(self.id) + "_" + smart_str(self.paper)

    def save(self, force_insert=False, force_update=False, using=None,
             update_fields=None):
        super(Datasets, self).save(force_insert, force_update, using, update_fields)
        self.complete_info(force_insert, force_update, using, update_fields)

    def complete_info(self, force_insert=False, force_update=False, using=None, update_fields=None):

        full_text = ""

        treatments = self.treats.all()
        if treatments:
            treatment_first = treatments[0]
            assert isinstance(treatment_first, Samples)
            treatment_first.update_a_dataset(self)

        for a_sample in self.treats.all():
            full_text += "\n".join([str(i) for i in [a_sample.unique_id, a_sample.name]])
        if self.factor:
            full_text += "\n" + self.factor.name
        if self.cell_line:
            full_text += "\n" + self.cell_line.name
        if self.cell_pop:
            full_text += "\n" + self.cell_pop.name
        if self.cell_type:
            full_text += "\n" + self.cell_type.name
        if self.strain:
            full_text += "\n" + self.strain.name

        if self.tissue_type:
            full_text += "\n" + self.tissue_type.name

        if self.disease_state:
            full_text += "\n" + self.disease_state.name

        if self.species:
            full_text += "\n" + self.species.name
        if self.paper:
            full_text += "\n".join(str(i) for i in [self.paper.authors, self.paper.title, self.paper.pmid])

        self.full_text = full_text

        if self.conts and self.conts.all() and self.conts.all()[0].unique_id == "clear controls":
            self.conts.clear()
        super(Datasets, self).save(force_insert, force_update, using, update_fields)



    journal_impact_factor.admin_order_field = 'paper__journal__impact_factor'

Datasets._meta._donotSerialize = ['user']
Datasets._meta._virtualfields = ['factor', 'platform', 'species',
                                 'assembly', 'cell_type', 'cell_line',
                                 'cell_pop', 'tissue_type', 'strain',
                                 'condition', 'disease_state',
                                 #not sure why treats and conts aren't in
                                 #_meta.fields, but pull them in too
                                 'treats', 'conts'
]

def samples_added(sender, **kwargs):
    if kwargs["action"].startswith("post"):
        kwargs["instance"].complete_info()


m2m_changed.connect(samples_added, sender=Datasets.treats.through)

class Samples(DCModel):
    """a table to store all of the sample information: a sample is a
    set of one or more datasets; it is associated with a paper (one paper to
    many samples)"""

    class Meta:
        verbose_name_plural = 'samples'

        #MIG_NOTE: can't rely on GSE!

    #def upload_factory(sub_dir):
    #    """a factory for generating upload_to_path fns--e.g. use to generate
    #    the various sub-directories we use to store the info associated w/
    #    a sample"""

    #    def upload_to_path(self, filename):
    #        """Returns the upload_to path for this sample.
    #        NOTE: we are going to store the files by the paper's gseid and
    #        the sample id, e.g. gseid = GSE20852, sample id = 578
    #        is going to be stored in: data/samples/gse20/852/578/[peak,wig,etc]
    #        """
    #        return os.path.join('data', 'samples',
    #                            'gse%s' % self.paper.gseid[3:5],
    #                            self.paper.gseid[5:], str(self.id), sub_dir,
    #                            filename)

    #    return upload_to_path

        #user = curator/uploader of this sample

    user = models.ForeignKey(User, null=True, blank=True, default=None, on_delete=models.CASCADE)
    paper = models.ForeignKey('Papers', null=True, blank=True, default=None, on_delete=models.SET_NULL)

    unique_id = models.CharField(max_length=255, null=True, blank=True, default="")
    #comma sep list of other identifiers
    other_ids = models.CharField(max_length=255, null=True, blank=True, default="")
    series_id = models.CharField(max_length=255, null=True, blank=True, default="")
    #Name comes from "title" in the geo sample information
    name = models.CharField(max_length=255, null=True, blank=True, default="")

    #RAW FILES assoc. w/ sample--i.e. FASTQ, and then when aligned --> BAM;
    #DELETE fastq when bam is generated
    #fastq_file = models.FileField(upload_to=upload_factory("fastq"),
    #                              null=True, blank=True, max_length=1024)
    fastq_file_url = models.CharField(max_length=255,
                                      null=True, blank=True)
    #fastq_file_size = models.CharField(max_length=255, null=True, blank=True, default="")
    #bam_file = models.FileField(upload_to=upload_factory("bam"),
    #                            null=True, blank=True, max_length=1024)

    #META information
    factor = models.ForeignKey('Factors', null=True, blank=True, default=None,on_delete=models.SET_NULL)
    platform = models.ForeignKey('Platforms',
                                 null=True, blank=True, default=None,on_delete=models.SET_NULL)
    species = models.ForeignKey('Species',
                                null=True, blank=True, default=None,on_delete=models.SET_NULL)
    assembly = models.ForeignKey('Assemblies',
                                 null=True, blank=True, default=None,on_delete=models.SET_NULL)
    #in description, we can add additional info e.g. protocols etc
    description = models.TextField(null=True, blank=True, default="")
    cell_type = models.ForeignKey('CellTypes',
                                  null=True, blank=True, default=None,on_delete=models.SET_NULL)
    cell_line = models.ForeignKey('CellLines',
                                  null=True, blank=True, default=None, on_delete=models.SET_NULL)
    cell_pop = models.ForeignKey('CellPops',
                                 null=True, blank=True, default=None,on_delete=models.SET_NULL)
    strain = models.ForeignKey('Strains',
                               null=True, blank=True, default=None,on_delete=models.SET_NULL)
    condition = models.ForeignKey('Conditions',
                                  null=True, blank=True, default=None,on_delete=models.SET_NULL)
    disease_state = models.ForeignKey('DiseaseStates',
                                      null=True, blank=True, default=None,on_delete=models.SET_NULL)
    tissue_type = models.ForeignKey('TissueTypes', null=True, blank=True,
                                    default=None,on_delete=models.SET_NULL)
    antibody = models.ForeignKey('Antibodies', null=True, blank=True,
                                 default=None, related_name='antibody',on_delete=models.SET_NULL)

    #curator = the person who double checks the info
    curator = models.ForeignKey(User, null=True, blank=True, default=None,
                                related_name="curator",on_delete=models.SET_NULL)

    status = models.CharField(max_length=255, choices=SAMPLE_STATUS,
                              null=True, blank=True, default="new")
    comments = models.TextField(null=True, blank=True, default="")

    dc_collect_date = models.DateTimeField(null=True, blank=True, default=None)
    dc_upload_date = models.DateTimeField(blank=True, null=True, default=None)

    geo_last_update_date = models.DateTimeField(blank=True, null=True, default=None)
    geo_release_date = models.DateTimeField(blank=True, null=True, default=None)

    re_check = models.TextField(null=True, blank=True, default=None)
    biological_source = models.CharField(max_length=255, null=True, blank=True, default=None)
    qc_judge = models.CharField(max_length=255, null=True, blank=True, default=None)

    is_correcting = models.BooleanField(default=False)

    def __str__(self):
        return smart_str(
            "_".join([self.unique_id, self.factor.name if self.factor else "", self.name if self.name else "", self.species.name if self.species else "NA"]))


    def save(self, *args, **kwargs):
        super(Samples, self).save(*args, **kwargs)
        self.update_related_datasets()

    def update_related_datasets(self):
        for a_dataset in self.TREATS.all():
            self.update_a_dataset(a_dataset)

    def update_a_dataset(self, dataset_obj):
        # Update in place
        dataset_obj.species = self.species
        dataset_obj.factor = self.factor
        dataset_obj.cell_line = self.cell_line
        dataset_obj.cell_type = self.cell_type
        dataset_obj.cell_pop = self.cell_pop
        dataset_obj.strain = self.strain
        dataset_obj.disease_state = self.disease_state
        dataset_obj.tissue_type = self.tissue_type
        dataset_obj.paper = self.paper
        super(Datasets, dataset_obj).save() 
    

Samples._meta._donotSerialize = ['user', 'curator']


class Platforms(DCModel):
    """Platforms are the chips/assemblies used to generate the dataset.
    For example, it can be an Affymetrix Human Genome U133 Plus 2.0 Array,
    i.e. GPLID = GPL570
    The fields are:
    name- name of the platform
    gplid - GEO Platform ID
    experiment type- Choice: ChIP-Chip/ChIP-Seq
    """
    gplid = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)
    name = models.CharField(max_length=255, null=True, blank=True, default="")
    technology = models.CharField(max_length=255, null=True, blank=True, default="")
    company = models.CharField(max_length=255, null=True, blank=True, default="")
    experiment_type = models.CharField(max_length=10, null=True, blank=True,
                                       default="",
                                       choices=EXPERIMENT_TYPE_CHOICES)

    def __str__(self):
        return smart_str(self.name)


class Aliases(DCModel):
    """Aliases of factors, e.g. p300, CBP, etc."""

    class Meta:
        verbose_name_plural = 'aliases'
    name = models.CharField(max_length=255)
    factor = models.ForeignKey('Factors', null=True, blank=True, default=None,related_name='aliases', on_delete=models.SET_NULL)

    def __str__(self):
        return smart_str(self.name)


class Factors(DCModel):
    """The factors applied to the sample, e.g. PolII, H3K36me3, etc."""

    class Meta:
        verbose_name_plural = 'factors'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    #antibody = models.CharField(max_length=255, blank=True)
    type = models.CharField(max_length=255, choices=FACTOR_TYPES,null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    def __str__(self):
        return smart_str(self.name)


class CellTypes(DCModel):
    """Sample's tissue/cell type, e.g. embryonic stem cell, b lymphocytes, etc.
    """

    class Meta:
        verbose_name_plural = 'cell types'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)

    #tissue_type = models.CharField(max_length=255, blank=True)
    def __str__(self):
        return smart_str(self.name)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)

class CellLines(DCModel):
    """Sample's cell lines.  I really don't know what distinguishes
    cell lines from cell populations or strains and mutations, but i'm going
    to create the tables just to be inclusive
    """

    class Meta:
        verbose_name_plural = 'cell lines'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)

    def __str__(self):
        return smart_str(self.name)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)


class CellPops(DCModel):
    class Meta:
        verbose_name_plural = 'cell populations'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)

    def __str__(self):
        return smart_str(self.name)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)

class Strains(DCModel):
    class Meta:
        verbose_name_plural = 'strains'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)

    def __str__(self):
        return smart_str(self.name)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)

class DiseaseStates(DCModel):
    """Information field for datasets"""

    class Meta:
        verbose_name_plural = 'disease states'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)
    def __str__(self):
        return smart_str(str(self.name))
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)



class Antibodies(DCModel):
    """Antibodies"""

    class Meta:
        verbose_name_plural = 'antibodies'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)
    def __str__(self):
        return smart_str(self.name)
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)


class TissueTypes(DCModel):
    """Tissue Types"""

    class Meta:
        verbose_name_plural = 'tissues'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    def __str__(self):
        return smart_str(self.name)
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)


class Conditions(DCModel):
    """Experiment/sample conditions, e.g. PTIP-knockout, wild-type"""

    class Meta:
        verbose_name_plural = 'conditions'
    comments = models.CharField(max_length=255, null=True, blank=True, default="")
    name = models.CharField(max_length=255)
    status = models.CharField(max_length=255, null=True, blank=True, default='new', choices=GENERAL_STATUS)

    def __str__(self):
        return smart_str(self.name)
    aliases =   models.CharField(max_length=255, null=True, blank=True, default=None)

class Journals(DCModel):
    """Journals that the papers are published in"""

    class Meta:
        verbose_name_plural = 'journals'

    name = models.CharField(max_length=255, null=True, blank=True, default="")
    issn = models.CharField(max_length=9, null=True, blank=True, default="")
    impact_factor = models.FloatField(default=0.0, null=True)

    def __str__(self):
        return smart_str(self.name)


class PaperSubmissions(DCModel):
    """Public paper submission page
    we collect the ip address of the submitter just in case of malicious usr
    pmid - pubmed id,
    gseid - GEO series id,
    NOTE: either pmid or gseid must be submitted
    status- see submission status,
    user- curator who last handled this submission
    ip_addr - the ip address of the submitter
    submitter_name - optional name of the submitter
    comments - any comments a currator might attach to the submission
    """
    #MIG_NOTE: need to change gseid to maybe unique_id or add something for
    #SRA/ENCODE etc
    pmid = models.IntegerField(default=0)
    gseid = models.CharField(max_length=8, blank=True)
    status = models.CharField(max_length=255, choices=SUBMISSION_STATUS)
    user = models.ForeignKey(User, null=True, blank=True, default=None, on_delete=models.CASCADE)
    ip_addr = models.CharField(max_length=15)
    submitter_name = models.CharField(max_length=255, blank=True)
    comments = models.TextField(blank=True)


PaperSubmissions._meta._donotSerialize = ['user']


class Species(DCModel):
    name = models.CharField(max_length=255, null=True, blank=True, default="")

    def __str__(self):
        return smart_str(self.name)


class Assemblies(DCModel):
    name = models.CharField(max_length=255, null=True, blank=True, default="")
    pub_date = models.DateField(blank=True)

    def __str__(self):
        return smart_str(self.name)


class UserProfiles(DCModel):
    """We want to add additional fields to the auth user model.  So creating
    this UserProfile model is the django way of doing it.
    ref: http://scottbarnham.com/blog/2008/08/21/extending-the-django-user-model-with-inheritance/
    NOTE: in the ref, the guy explains how to do it through model inheritance,
    but get_profile is now a django idiom that i decided to use it instead.
    """
    user = models.OneToOneField(User, unique=True, related_name='profile', on_delete=models.CASCADE)
    #which team is the user on, paper team or data team
    team = models.CharField(max_length=255, choices=TEAMS, blank=True, 
                            null=True, default=None)


class CistromeUser(models.Model):
    id = models.IntegerField(primary_key=True)
    create_time = models.DateTimeField()
    update_time = models.DateTimeField()
    email = models.CharField(max_length=255)
    password = models.CharField(max_length=255)
    external = models.IntegerField(null=True, blank=True)
    deleted = models.IntegerField(null=True, blank=True)
    purged = models.IntegerField(null=True, blank=True)
    username = models.CharField(max_length=255, blank=True)
    form_values_id = models.IntegerField(null=True, blank=True)
    disk_usage = models.DecimalField(null=True, max_digits=16, decimal_places=0, blank=True)
    active = models.IntegerField(null=True, blank=True)
    activation_token = models.CharField(max_length=64, blank=True)
    class Meta:
        db_table = 'galaxy_user'
        in_db = "cistromeap"


class Corrections(models.Model): #Xin_text
    id = models.AutoField(primary_key=True)
    dc_id = models.ForeignKey('Samples', to_field='id', null=True, blank=True, on_delete=models.CASCADE)
    treatment = models.CharField(max_length=255, null=True, blank=True, default=None)
    species = models.CharField(max_length=255, null=True, blank=True, default=None)
    pmid = models.CharField(max_length=255, null=True, blank=True, default=None)
    factor = models.CharField(max_length=255, null=True, blank=True, default=None)
    factor_type = models.CharField(max_length=255, null=True, blank=True, default=None)
    cell_line = models.CharField(max_length=255, null=True, blank=True, default=None)
    cell_type = models.CharField(max_length=255, null=True, blank=True, default=None)
    cell_pop = models.CharField(max_length=255, null=True, blank=True, default=None)
    strain = models.CharField(max_length=255, null=True, blank=True, default=None)
    tissue = models.CharField(max_length=255, null=True, blank=True, default=None)
    disease = models.CharField(max_length=255, null=True, blank=True, default=None)
    comment = models.CharField(max_length=255, null=True, blank=True, default=None)
    email = models.CharField(max_length=255, null=True, blank=True, default=None)
    time = models.DateTimeField(auto_now_add=True)
    last_treatment = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_species = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_pmid = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_factor = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_factor_type = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_cell_line = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_cell_type = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_cell_pop = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_strain = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_tissue = models.CharField(max_length=255, null=True, blank=True, default=None)
    last_disease = models.CharField(max_length=255, null=True, blank=True, default=None)
    status = models.IntegerField(blank=False, default=0)  # 0 means unchecked, 1 means checked, -1 means deleted, 2 means destroyed

    class Meta:
        verbose_name_plural = 'corrections'



# UserProfiles._meta._donotSerialize = ['user']

# class SampleDhsStats(DCModel):
#     """Stats about the sample"""
#     sample = models.ForeignKey('Samples', unique=True)
#     total_peaks = models.IntegerField(default=0)
#     peaks_in_dhs = models.IntegerField(default=0)
