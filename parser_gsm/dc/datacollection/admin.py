import re
from adminactions.mass_update import MassUpdateForm
from django.contrib import admin

from django import forms
from django.http import HttpResponse
from django_select2.widgets import Select2Widget, HeavySelect2MultipleWidget
from django.utils.translation import ugettext_lazy as _

from datacollection.admin_utils import SamplesSelect2View
from models import *

# from datacollection.widgets import *


import adminactions.actions as actions

from django.contrib.admin import site
# register all adminactions
actions.add_to_site(site)



class ImpactFactorFilter(admin.SimpleListFilter):
    # Human-readable title which will be displayed in the
    # right admin sidebar just above the filter options.
    title = _('impact factor')

    # Parameter for the filter that will be used in the URL query.
    parameter_name = 'impact'

    def lookups(self, request, model_admin):
        """
        Returns a list of tuples. The first element in each
        tuple is the coded value for the option that will
        appear in the URL query. The second element is the
        human-readable name for the option that will appear
        in the right sidebar.
        """
        return (
            ('>5', _('larger than 5')),
            ('>10', _('larger than 10')),
        )

    def queryset(self, request, queryset):
        """
        Returns the filtered queryset based on the value
        provided in the query string and retrievable via
        `self.value()`.
        """
        # Compare the requested value (either '80s' or '90s')
        # to decide how to filter the queryset.
        if self.value() == '>5':
            return queryset.filter(paper__journal__impact_factor__gt=5)
        if self.value() == '>10':
            return queryset.filter(paper__journal__impact_factor__gt=10)



class EmptyFilter(admin.SimpleListFilter):
    # Human-readable title which will be displayed in the
    # right admin sidebar just above the filter options.
    title = _('empty')

    # Parameter for the filter that will be used in the URL query.
    parameter_name = 'empty'

    def lookups(self, request, model_admin):
        """
        Returns a list of tuples. The first element in each
        tuple is the coded value for the option that will
        appear in the URL query. The second element is the
        human-readable name for the option that will appear
        in the right sidebar.
        """
        return (
            ('empty', _('no treatments')),
        )

    def queryset(self, request, queryset):
        """
        Returns the filtered queryset based on the value
        provided in the query string and retrievable via
        `self.value()`.
        """
        # Compare the requested value (either '80s' or '90s')
        # to decide how to filter the queryset.
        if self.value() == 'empty':
            return queryset.filter(treats__isnull=True)






class DatasetInline(admin.TabularInline):
    model = Datasets
    fields = ['user', 'paper', 'treats', 'conts', 'status']
    raw_id_fields = ['treats', 'conts']


class SampleInline(admin.TabularInline):
    model = Samples
    fields = ['paper', 'unique_id', 'series_id', 'name']


class TreatInline(admin.TabularInline):
    formfield_overrides = {
        models.ForeignKey: {'widget': Select2Widget(attrs={"style": "width: 300px;"})},
    }
    verbose_name = "treatment samples"
    name = verbose_name
    model = Datasets.treats.through
    extra = 1

    def formfield_for_foreignkey(self, db_field, request=None, **kwargs):
        if request._dataset_:
            series_id = request._dataset_.treats.all()[0].series_id
        else:
            series_id = None
        field = super(TreatInline, self).formfield_for_foreignkey(db_field, request, **kwargs)

        if not series_id:
            return field

        if db_field.name == 'samples':
            field.queryset = Samples.objects.filter(series_id=series_id)
        elif db_field.name == 'datasets':
            field.queryset = Datasets.objects.filter(treats__series_id=series_id)
        else:
            dn = db_field.name
            raise
        return field


class ContInline(admin.TabularInline):
    formfield_overrides = {
        models.ForeignKey: {'widget': Select2Widget(attrs={"style": "width: 300px;"})},
    }
    verbose_name = "control samples"
    name = verbose_name
    model = Datasets.conts.through
    extra = 1

    def formfield_for_foreignkey(self, db_field, request=None, **kwargs):
        if request._dataset_:
            series_id = request._dataset_.treats.all()[0].series_id
        else:
            series_id = None

        field = super(ContInline, self).formfield_for_foreignkey(db_field, request, **kwargs)

        if not series_id:
            return field

        if db_field.name == 'samples':
            field.queryset = Samples.objects.filter(series_id=series_id)
        elif db_field.name == 'datasets':
            field.queryset = Datasets.objects.filter(conts__series_id=series_id)
        else:
            dn = db_field.name
            raise
        return field


class PaperAdmin(admin.ModelAdmin):
    list_display = ['pmid', 'title', 'journal', 'date_collected', 'pub_date', 'status', 'reference', 'authors']
    list_filter = ['journal', 'date_collected', 'pub_date', 'status']
    search_fields = ['title', 'abstract', 'journal__name']
    list_per_page = 50
    inlines = [DatasetInline]
    list_max_show_all = 5000

    def suit_row_attributes(self, obj):
        css_class = {
            'complete': 'success',
            'error': 'error',
            'transfer': 'info',
            'datasets': 'warning',
            'imported': 'warning'
        }.get(obj.status)
        if css_class:
            return {'class': css_class}


class DatasetsHeavySelect2MultipleWidget(HeavySelect2MultipleWidget):
    def __init__(self, **kwargs):
        return super(DatasetsHeavySelect2MultipleWidget, self).__init__(data_view=SamplesSelect2View,
                                                                        data_url="http://cistrome.org/dc/select2sample/",
                                                                        **kwargs)


class DatasetMassUpdateForm(MassUpdateForm):
    conts = forms.ModelMultipleChoiceField(queryset=Samples.objects.all(),
                                           widget=DatasetsHeavySelect2MultipleWidget, required=False)

    class Meta:
        model = Datasets

    def model_fields(self):
        return [field for field in super(DatasetMassUpdateForm, self).model_fields() if
                field.name not in ['treats'] and not field.name.endswith(
                    "file") and not field.name.startswith("rep")]


class DatasetAddForm(forms.ModelForm):
    treats = forms.ModelMultipleChoiceField(queryset=Samples.objects.all(), widget=DatasetsHeavySelect2MultipleWidget)
    conts = forms.ModelMultipleChoiceField(queryset=Samples.objects.all(),
                                           widget=DatasetsHeavySelect2MultipleWidget, required=False)

    class Meta:
        model = Datasets


class DatasetAdmin(admin.ModelAdmin):
    class Media:
        js = ("new_admin.js",)
    mass_update_form = DatasetMassUpdateForm

    formfield_overrides = {
        models.ForeignKey: {'widget': Select2Widget(attrs={"style": "width: 300px;"})},
    }
    list_display = ['custom_id', 'paper', 'status', 'journal_name', 'journal_impact_factor', 'factor', 'custom_treats',
                    'custom_conts', 'custom_gse', 'action']
    list_filter = ['status', 'paper__journal', 'treats__factor__name', 'treats__factor__type', ImpactFactorFilter, EmptyFilter]
    search_fields = ['id', 'paper__title', 'paper__pmid', 'paper__journal__name', 'treats__factor__name',
                     'treats__unique_id',
                     'conts__unique_id', 'treats__series_id', 'conts__series_id']
    list_per_page = 100
    list_display_links = ['action']
    list_max_show_all = 5000



    def change_view(self, request, object_id, form_url='', extra_context=None):
        self.inlines = [TreatInline, ContInline, ]
        self.fields = ['user', 'paper', 'date_created', 'status', 'comments','full_text' ]
        return super(DatasetAdmin, self).change_view(request, object_id, form_url, extra_context)

    def add_view(self, request, form_url='', extra_context=None):
        self.inlines = []
        self.form = DatasetAddForm
        self.fields = ['user', 'paper', 'treats', 'conts', 'date_created', 'status', 'comments', ]
        return super(DatasetAdmin, self).add_view(request, form_url, extra_context)


    def suit_row_attributes(self, obj):
        css_class = {
            'validated': 'success',
            'complete':'success'
        }.get(obj.status)
        if css_class:
            return {'class': css_class}

    def get_form(self, request, obj=None, **kwargs):
    # just save sample reference for future processing in Treatment & Control Inline
        request._dataset_ = obj
        return super(DatasetAdmin, self).get_form(request, obj, **kwargs)


    def action(self, obj):
        return "Change"

    def custom_id(self, obj):
        return '<a href="http://cistrome.org/finder/util/meta?did=%s" target="_blank">%s</a>' % (obj.id, obj.id)


    geo_link_template = '<a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">%s</a>' \
                        '<a href="http://cistrome.org/dc/new_admin/datacollection/samples/?q=%s" target="_blank"><i class="icon-search"></i></a></br> %s'

    def custom_treats(self, obj):
        treats_list = obj.treats.all().order_by('unique_id')
        return "</br>".join(
            [DatasetAdmin.geo_link_template % (i.unique_id,
                                               i.unique_id,
                                               " ".join(re.findall(r"\d+", i.unique_id)),
                                               i.name) for i in
             treats_list])

    gse_link_template = '<a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE%s" target="_blank">GSE%s</a>' \
                        '<a href="http://cistrome.org/dc/new_admin/datacollection/samples/?q=%s" target="_blank"><i class="icon-search"></i></a></br>' \
                        '<a href="http://cistrome.org/dc/new_admin/datacollection/datasets/?q=GSE%s" target="_blank"><i class="icon-magnet"></i></a></br>'

    def custom_gse(self, obj):
        treats_list = obj.treats.all()
        if treats_list:
            try:
                gse_id = json.loads(treats_list[0].other_ids)['gse'].strip()[-5:]
                return DatasetAdmin.gse_link_template % (gse_id, gse_id, gse_id, gse_id)
            except:
                return ""
        return ""

    def custom_conts(self, obj):
        conts_list = obj.conts.all().order_by('unique_id')
        return "</br>".join(
            [DatasetAdmin.geo_link_template % (i.unique_id,
                                               i.unique_id,
                                               " ".join(re.findall(r"\d+", i.unique_id)),
                                               i.name) for i in
             conts_list])


    custom_id.allow_tags = True
    custom_id.short_description = 'ID'
    custom_id.ordering = 'id'

    custom_treats.allow_tags = True
    custom_treats.short_description = "Treatment IDs"

    custom_conts.allow_tags = True
    custom_conts.short_description = "Control IDs"

    custom_gse.allow_tags = True
    custom_gse.short_description = "GSE ID"


class SampleMassUpdateForm(MassUpdateForm):

    class Meta:

        model = Samples
        exclude = ("dc_collect_date", "dc_upload_date","geo_last_update_date","geo_release_date","user","paper","fastq_file","fastq_file_url","bam_file","platform","paper")
        widgets = {
            # 'paper': Select2Widget(attrs={"style": "width: 300px;"}),
            'factor': Select2Widget(attrs={"style": "width: 300px;"}),
            'cell_line': Select2Widget(attrs={"style": "width: 300px;"}),
            'cell_type': Select2Widget(attrs={"style": "width: 300px;"}),
            'cell_pop': Select2Widget(attrs={"style": "width: 300px;"}),
            'condition': Select2Widget(attrs={"style": "width: 300px;"}),
            'strain': Select2Widget(attrs={"style": "width: 300px;"}),
            'disease_state': Select2Widget(attrs={"style": "width: 300px;"}),
            'tissue_type': Select2Widget(attrs={"style": "width: 300px;"}),
            'antibody': Select2Widget(attrs={"style": "width: 300px;"}),
            # 'platform': Select2Widget(attrs={"style": "width: 300px;"}),
            # 'species': Select2Widget(attrs={"style": "width: 300px;"}),
            # 'assembly': Select2Widget(attrs={"style": "width: 300px;"}),
        }

class CellInfoFilter(admin.SimpleListFilter):
    title = _("cell info")
    parameter_name = 'cellinfo'
    def lookups(self, request, model_admin):
        return (
            ('lack', _('lacks cell information')),
            ('ok', _('has cell information'))
        )
    def queryset(self, request, queryset):
        if self.value() == 'lack':
            return queryset.filter(cell_line__name=None, cell_type__name=None, cell_pop__name=None, tissue_type__name=None)
        if self.value() == 'ok':
            return queryset.exclude(cell_line__name=None, cell_type__name=None, cell_pop__name=None, tissue_type__name=None)

class FactorInfoFilter(admin.SimpleListFilter):
    title = _("factor info")
    parameter_name = 'factorinfo'
    def lookups(self, request, model_admin):
        return (
            ('lack', _('lacks factor information')),
            ('ok', _('has cell information'))
        )
    def queryset(self, request, queryset):
        if self.value() == 'lack':
            return queryset.filter(factor__name=None)
        if self.value() == 'ok':
            return queryset.exclude(factor__name=None)


class GroupingInfoFilter(admin.SimpleListFilter):
    title = _("grouping info")
    parameter_name = 'groupinginfo'
    def lookups(self, request, model_admin):
        return (
            ('treatment', _('grouped as treatment')),
            ('control', _('grouped as control')),
            ('orphan', _('orphan (no grouping info)'))
        )

    def queryset(self, request, queryset):
        if self.value() == 'treatment':
            return queryset.exclude(TREATS__isnull=True)
        if self.value() == 'control':
            return queryset.exclude(CONTS__isnull=True)
        if self.value() == 'orphan':
            return queryset.filter(TREATS__isnull=True, CONTS__isnull=True)


class RecheckInfoFilter(admin.SimpleListFilter):
    title = _("recheck info")
    parameter_name = 'recheckinfo'
    def lookups(self, request, model_admin):
        return (
            ('need', _('need recheck')),
        )
    def queryset(self, request, queryset):
        if self.value() == 'need':
            return queryset.exclude(re_check=None).exclude(re_check="{}").exclude(re_check="")



class SampleAdmin(admin.ModelAdmin):
    class Media:
        js = ("new_admin.js",)
    mass_update_form = SampleMassUpdateForm

    formfield_overrides = {
        models.ForeignKey: {'widget': Select2Widget(attrs={"style": "width: 300px;"})},
    }
    list_display = ['id', 'unique_id_url', 'other_id', 'status', 'custom_name', 'factor', 'species', 'cell_category',
                    'cell_source',
                    'condition',
                    'custom_antibody', 'custom_description', 'paper','action', 'custom_recheck']
    search_fields = ['id', 'unique_id', 'other_ids', 'factor__name', 'species__name', 'cell_type__name',
                     'cell_line__name',
                     'cell_pop__name', 'strain__name', 'name',
                     'condition__name',
                     'disease_state__name', 'tissue_type__name', 'antibody__name', 'description', 'series_id']
    list_display_links = ['action']
    list_filter = ['status', 'species__name', 'factor__name', 'factor__type',  CellInfoFilter, FactorInfoFilter, GroupingInfoFilter, RecheckInfoFilter ]
    list_per_page = 200


    def suit_row_attributes(self, obj, request):
        css_class = {
            'imported': 'success',
            'checked': 'success',
            'new': 'warning',
            "ignored": 'info',
        }.get(obj.status)
        if css_class:
            return {'class': css_class}

    def custom_recheck(self, obj):
        ret = ""
        if obj.re_check:
            recheck_dict = json.loads(obj.re_check)
            for k in recheck_dict:
                ret += "<strong>%s</strong>: %s<br>" % (k.title(), recheck_dict[k])
        return ret

    def custom_name(self, obj):
        if obj.name:
            return obj.name.replace("_", " ")

    def unique_id_url(self, obj):

        return '<a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s" target="_blank">%s</a>' % (
            obj.unique_id, obj.unique_id)


    def other_id(self, obj):
        ret = ""
        try:
            other_ids = json.loads(obj.other_ids)
        except:
            return ret
        if other_ids['pmid']:
            pmid = other_ids['pmid'].strip()
            ret += '<a href="http://www.ncbi.nlm.nih.gov/pubmed/%s" target="_blank">P%s</a><br>' % (pmid, pmid)
        if other_ids['gse']:
            gse = other_ids['gse'].strip()[-5:]
            ret += '<br><a href="http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE%s" target="_blank">GSE%s</a>' % (
                gse, gse)
        return ret

    def action(self, obj):
        return "Change"

    def custom_antibody(self, obj):
        return obj.antibody

    def custom_description(self, obj):
        ret = ""
        if obj.description:
            desc_dict = json.loads(obj.description)
            for k in desc_dict:
                ret += "<strong>%s</strong>: %s<br>" % (k.title(), desc_dict[k])
        return ret

    def cell_source(self, obj):
        ret = ""
        if obj.cell_line:
            ret += "<strong>Cell Line</strong>: %s<br>" % obj.cell_line
        if obj.cell_pop:
            ret += "<strong>Cell Pop</strong>: %s<br>" % obj.cell_pop
        if obj.strain:
            ret += "<strong>Strain</strong>: %s<br>" % obj.strain
        return ret

    def cell_category(self, obj):
        ret = ""
        if obj.cell_type:
            ret += "<strong>Cell Type</strong>: %s<br>" % obj.cell_type

        if obj.disease_state:
            ret += "<strong>Disease</strong>: %s<br>" % obj.disease_state

        if obj.tissue_type:
            ret += "<strong>Tissue</strong>: %s<br>" % obj.tissue_type

        return ret

    custom_antibody.allow_tags = True
    custom_antibody.short_description = 'Antibody'
    custom_antibody.admin_order_field = 'antibody'
    custom_description.allow_tags = True
    custom_description.short_description = "Description"
    custom_recheck.allow_tags = True
    custom_recheck.short_description = "Recheck"
    custom_name.short_description = "Name"
    custom_name.admin_order_field = 'name'
    cell_source.allow_tags = True
    cell_category.allow_tags = True

    unique_id_url.short_description = 'Sample ID'
    unique_id_url.allow_tags = True
    unique_id_url.admin_order_field = 'unique_id'
    other_id.short_description = 'Other ID'
    other_id.allow_tags = True

    ordering = ['-id']
    list_max_show_all = 5000


class JournalAdmin(admin.ModelAdmin):
    list_display = ['name', 'issn', 'impact_factor']
    list_per_page = 500


class FactorAdmin(admin.ModelAdmin):
    list_display = ['id', 'name', 'custom_aliases', 'custom_type', 'status', ]
    list_filter = ['type','status']
    inlines = [SampleInline, ]
    list_per_page = 100

    def suit_row_attributes(self, obj):
        css_class = {
            'new': 'warning'
        }.get(obj.status)
        if css_class:
            return {'class': css_class}


    def custom_type(self, obj):
        defined_choice = obj.get_type_display()
        if defined_choice:
            return defined_choice
        else:
            return obj.type

    def custom_aliases(self, obj):
        aliases_list = []
        for i in obj.aliases.all():
            aliases_list.append(i.name)
        return ", ".join(aliases_list)

    custom_type.short_description = 'Type'
    custom_type.admin_order_field = "type"

    search_fields = ['id', 'name']

    list_max_show_all = 5000


class CellInfoAdmin(admin.ModelAdmin):
    # inlines = [SampleInline, ]
    list_per_page = 500
    list_display = ['id', 'name', 'status','aliases']
    search_fields = ['id', 'name', 'aliases']
    list_max_show_all = 5000
    list_filter = ['status']
    def suit_row_attributes(self, obj):
        css_class = {
            'new': 'warning'
        }.get(obj.status)
        if css_class:
            return {'class': css_class}

class TissueAdmin(CellInfoAdmin):
    pass

class CelllineAdmin(CellInfoAdmin):
    pass


class CelltypeAdmin(CellInfoAdmin):
    pass


class CellpopAdmin(CellInfoAdmin):
    pass


class DiseaseAdmin(CellInfoAdmin):
    pass


class StrainAdmin(CellInfoAdmin):
    pass


class AntibodyAdmin(CellInfoAdmin):
    pass


class ConditionAdmin(CellInfoAdmin):
    pass


class AliasAdmin(admin.ModelAdmin):
    formfield_overrides = {
        models.ForeignKey: {'widget': Select2Widget(attrs={"style": "width: 300px;"})},
        }
    list_display = ['id', 'name', 'factor']
    search_fields = ['id', 'name']
    list_per_page = 100
    search_field = ['name', 'factor__name']
    list_max_show_all = 5000


admin.site.register(Aliases, AliasAdmin)
admin.site.register(Papers, PaperAdmin)
admin.site.register(Samples, SampleAdmin)
admin.site.register(Datasets, DatasetAdmin)
admin.site.register(Journals, JournalAdmin)
admin.site.register(Factors, FactorAdmin)
admin.site.register(CellTypes, CelltypeAdmin)
admin.site.register(CellLines, CelllineAdmin)
admin.site.register(CellPops, CellpopAdmin)
admin.site.register(DiseaseStates, DiseaseAdmin)
admin.site.register(Strains, StrainAdmin)
admin.site.register(Antibodies, AntibodyAdmin)
admin.site.register(Conditions, ConditionAdmin)
admin.site.register(TissueTypes, TissueAdmin)

