# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'zi_control.ui'
#
# Created: Wed Mar 23 12:47:45 2016
#      by: PyQt4 UI code generator 4.10.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1138, 859)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.list_records = QtGui.QListView(self.centralwidget)
        self.list_records.setGeometry(QtCore.QRect(740, 290, 391, 141))
        self.list_records.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.list_records.setModelColumn(0)
        self.list_records.setObjectName(_fromUtf8("list_records"))
        self.layoutWidget = QtGui.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(730, 500, 411, 25))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.layoutWidget)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.txt_tag = QtGui.QLineEdit(self.layoutWidget)
        self.txt_tag.setObjectName(_fromUtf8("txt_tag"))
        self.horizontalLayout_2.addWidget(self.txt_tag)
        self.btn_clear_record = QtGui.QPushButton(self.layoutWidget)
        self.btn_clear_record.setObjectName(_fromUtf8("btn_clear_record"))
        self.horizontalLayout_2.addWidget(self.btn_clear_record)
        self.btn_save_record_to_disk = QtGui.QPushButton(self.layoutWidget)
        self.btn_save_record_to_disk.setObjectName(_fromUtf8("btn_save_record_to_disk"))
        self.horizontalLayout_2.addWidget(self.btn_save_record_to_disk)
        self.tabWidget = QtGui.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 340, 721, 471))
        self.tabWidget.setTabShape(QtGui.QTabWidget.Rounded)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab_scripts = QtGui.QWidget()
        self.tab_scripts.setObjectName(_fromUtf8("tab_scripts"))
        self.tree_scripts = QtGui.QTreeWidget(self.tab_scripts)
        self.tree_scripts.setEnabled(True)
        self.tree_scripts.setGeometry(QtCore.QRect(0, 10, 711, 291))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tree_scripts.sizePolicy().hasHeightForWidth())
        self.tree_scripts.setSizePolicy(sizePolicy)
        self.tree_scripts.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        self.tree_scripts.setHeaderHidden(False)
        self.tree_scripts.setObjectName(_fromUtf8("tree_scripts"))
        self.tree_scripts.header().setDefaultSectionSize(150)
        self.tree_scripts.header().setHighlightSections(True)
        self.progressBar = QtGui.QProgressBar(self.tab_scripts)
        self.progressBar.setGeometry(QtCore.QRect(230, 410, 471, 23))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName(_fromUtf8("progressBar"))
        self.btn_stop_script = QtGui.QPushButton(self.tab_scripts)
        self.btn_stop_script.setGeometry(QtCore.QRect(120, 410, 101, 23))
        self.btn_stop_script.setObjectName(_fromUtf8("btn_stop_script"))
        self.btn_start_script = QtGui.QPushButton(self.tab_scripts)
        self.btn_start_script.setGeometry(QtCore.QRect(20, 410, 101, 23))
        self.btn_start_script.setObjectName(_fromUtf8("btn_start_script"))
        self.list_scripts = QtGui.QListView(self.tab_scripts)
        self.list_scripts.setGeometry(QtCore.QRect(0, 320, 711, 61))
        self.list_scripts.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)
        self.list_scripts.setModelColumn(0)
        self.list_scripts.setObjectName(_fromUtf8("list_scripts"))
        self.tabWidget.addTab(self.tab_scripts, _fromUtf8(""))
        self.tab_monitor = QtGui.QWidget()
        self.tab_monitor.setObjectName(_fromUtf8("tab_monitor"))
        self.tree_monitor = QtGui.QTreeWidget(self.tab_monitor)
        self.tree_monitor.setEnabled(True)
        self.tree_monitor.setGeometry(QtCore.QRect(0, 10, 711, 431))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tree_monitor.sizePolicy().hasHeightForWidth())
        self.tree_monitor.setSizePolicy(sizePolicy)
        self.tree_monitor.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        self.tree_monitor.setHeaderHidden(False)
        self.tree_monitor.setObjectName(_fromUtf8("tree_monitor"))
        self.tree_monitor.header().setDefaultSectionSize(150)
        self.tree_monitor.header().setHighlightSections(True)
        self.tabWidget.addTab(self.tab_monitor, _fromUtf8(""))
        self.tab_settings = QtGui.QWidget()
        self.tab_settings.setObjectName(_fromUtf8("tab_settings"))
        self.tree_settings = QtGui.QTreeWidget(self.tab_settings)
        self.tree_settings.setEnabled(True)
        self.tree_settings.setGeometry(QtCore.QRect(0, 10, 711, 431))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tree_settings.sizePolicy().hasHeightForWidth())
        self.tree_settings.setSizePolicy(sizePolicy)
        self.tree_settings.setEditTriggers(QtGui.QAbstractItemView.AllEditTriggers)
        self.tree_settings.setHeaderHidden(False)
        self.tree_settings.setObjectName(_fromUtf8("tree_settings"))
        self.tree_settings.header().setDefaultSectionSize(150)
        self.tree_settings.header().setHighlightSections(True)
        self.tabWidget.addTab(self.tab_settings, _fromUtf8(""))
        self.tab_history = QtGui.QWidget()
        self.tab_history.setObjectName(_fromUtf8("tab_history"))
        self.tabWidget.addTab(self.tab_history, _fromUtf8(""))
        self.matplotlibwidget = MatplotlibWidget(self.centralwidget)
        self.matplotlibwidget.setGeometry(QtCore.QRect(30, 20, 400, 300))
        self.matplotlibwidget.setObjectName(_fromUtf8("matplotlibwidget"))
        self.list_history = QtGui.QListView(self.centralwidget)
        self.list_history.setGeometry(QtCore.QRect(470, 10, 631, 211))
        self.list_history.setObjectName(_fromUtf8("list_history"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1138, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(2)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow", None))
        self.txt_tag.setText(_translate("MainWindow", "some tag", None))
        self.btn_clear_record.setText(_translate("MainWindow", "Delete", None))
        self.btn_save_record_to_disk.setText(_translate("MainWindow", "Save to disk", None))
        self.tree_scripts.headerItem().setText(0, _translate("MainWindow", "Parameter/Script", None))
        self.tree_scripts.headerItem().setText(1, _translate("MainWindow", "Value", None))
        self.btn_stop_script.setText(_translate("MainWindow", "stop", None))
        self.btn_start_script.setText(_translate("MainWindow", "start", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_scripts), _translate("MainWindow", "Scripts", None))
        self.tree_monitor.headerItem().setText(0, _translate("MainWindow", "Parameter", None))
        self.tree_monitor.headerItem().setText(1, _translate("MainWindow", "Value", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_monitor), _translate("MainWindow", "Monitor", None))
        self.tree_settings.headerItem().setText(0, _translate("MainWindow", "Instrument", None))
        self.tree_settings.headerItem().setText(1, _translate("MainWindow", "Value", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_settings), _translate("MainWindow", "Settings", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_history), _translate("MainWindow", "History", None))

from matplotlibwidget import MatplotlibWidget
